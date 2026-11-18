#include "equilibrium/spdata_output.hpp"
#include "equilibrium/spdata_validation.hpp"
#include "spdata_fixture.hpp"
#include "test_support.hpp"

#include <array>
#include <fstream>
#include <iterator>
#include <limits>
#include <locale>
#include <sstream>
#include <string>

#if defined(__linux__) && !defined(__EMSCRIPTEN__)
#include <csignal>

#include <sys/resource.h>
#include <sys/wait.h>
#endif

namespace {
using equilibrium::test::check;

std::string contents(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(file), {}};
}

struct CommaDecimal : std::numpunct<char> {
    char do_decimal_point() const override { return ','; }
};

void check_file(std::size_t lst) {
    const auto data = spdata_fixture(lst);
    const auto path = "spdata_" + std::to_string(lst) + ".dat";
    const auto status = equilibrium::write_spdata(path, data);
    check(status.has_value(), "spdata writer succeeds");
    if (!status) return;
    std::istringstream input(contents(path));
    input.imbue(std::locale::classic());
    std::string line;
    std::getline(input, line);
    check(line == data.meta_data, "metadata is the first line");
    std::getline(input, line);
    check(line.size() == 16 && line.substr(8) == "   4   8",
          "I4 header matches MCT");
    std::istringstream dimensions(line);
    std::size_t nr = 0, nt = 0;
    dimensions >> nr >> nt;
    check(nr == data.lsp && nt == lst, "dimensions parse as separate tokens");

    const auto row = [&](std::span<const double> expected) {
        if (!std::getline(input, line) || line.size() != expected.size() * 18)
            throw std::runtime_error("incorrect spdata row length");
        for (std::size_t k = 0; k < expected.size(); ++k) {
            const auto cell = line.substr(18 * k, 18);
            if (cell.find('E') == std::string::npos ||
                cell.find(',') != std::string::npos ||
                std::stod(cell) != expected[k])
                throw std::runtime_error("spdata coefficient changed or moved");
        }
    };
    row(std::array{data.psi_wall, data.psi_sep});
    const std::array fields{&data.b, &data.r, &data.z, &data.j};
    const std::array profiles{&data.q, &data.f,       &data.i,
                              &data.p, &data.r_minor, &data.psi_t};
    for (std::size_t r = 0; r < nr; ++r) {
        for (const auto* field : fields) {
            for (std::size_t c = 0; c < 9; ++c) {
                for (std::size_t t = 0; t <= nt;) {
                    std::array<double, 4> values{};
                    std::size_t count = 0;
                    while (count < 4 && t <= nt) {
                        values[count++] = (*field)[(r * 9 + c) * nt + t % nt];
                        ++t;
                    }
                    row(std::span(values).first(count));
                }
            }
        }
        for (const auto* profile : profiles)
            row(std::span(*profile).subspan(r * 3, 3));
    }
    std::getline(input, line);
    check(line == "   0   0", "axisymmetric ripple header");
    row(std::array{2.0, 0.0, 0.0});
    row(std::array{0.0, 0.0});
    check(!std::getline(input, line),
          "complete spdata footer and no extra rows");
}

void failures() {
    const std::string path = "spdata_preserved.dat";
    { std::ofstream(path) << "original"; }
    const auto valid = spdata_fixture(5);
    const auto rejects = [&](const auto& change) {
        auto bad = valid;
        change(bad);
        check(!equilibrium::write_spdata(path, bad), "invalid record rejected");
        check(contents(path) == "original", "invalid output preserves target");
    };
    rejects([](auto& d) { d.lsp = 1; });
    rejects([](auto& d) { d.lsp = 10000; });
    rejects([](auto& d) { d.lst = 0; });
    rejects([](auto& d) { d.lst = 1000; });
    rejects([](auto& d) { d.lsp = std::numeric_limits<std::size_t>::max(); });
    rejects([](auto& d) { d.psi_wall = 0.0; });
    rejects(
        [](auto& d) { d.psi_sep = std::numeric_limits<double>::infinity(); });
    rejects([](auto& d) { d.meta_data += '\n'; });
    rejects([](auto& d) { d.meta_data += '\r'; });
    rejects([](auto& d) { d.meta_data += '\0'; });
    rejects([](auto& d) { d.r.pop_back(); });
    rejects([](auto& d) { d.psi_t.pop_back(); });
    rejects(
        [](auto& d) { d.j.back() = std::numeric_limits<double>::quiet_NaN(); });
    rejects(
        [](auto& d) { d.p.back() = std::numeric_limits<double>::infinity(); });
    rejects([](auto& d) { d.b[1] += 1.0; });
    rejects([](auto& d) { d.r[d.field_index(0, 3, 0)] = 1.0; });
    rejects([](auto& d) { std::fill_n(d.r.begin(), d.lst, 0.0); });
    check(!equilibrium::write_spdata("", valid), "empty path rejected");
    check(!equilibrium::write_spdata(path + std::string("\0tail", 5), valid),
          "embedded NUL in path rejected");
    check(!equilibrium::write_spdata("absent_directory/spdata.dat", valid),
          "open failure reported");
    std::filesystem::create_directory("spdata_target_directory");
    check(!equilibrium::write_spdata("spdata_target_directory", valid),
          "rename failure reported");
    std::filesystem::remove("spdata_target_directory");

#if defined(__linux__) && !defined(__EMSCRIPTEN__)
    const auto child = fork();
    if (child == 0) {
        std::signal(SIGXFSZ, SIG_IGN);
        const rlimit limit{1024, 1024};
        if (setrlimit(RLIMIT_FSIZE, &limit) != 0) _exit(2);
        _exit(equilibrium::write_spdata(path, valid) ? 3 : 0);
    }
    int status = 0;
    check(child > 0 && waitpid(child, &status, 0) == child &&
              WIFEXITED(status) && WEXITSTATUS(status) == 0,
          "storage write failure reported");
    check(contents(path) == "original", "storage failure preserves target");
#endif
    bool stray = false;
    for (const auto& entry : std::filesystem::directory_iterator("."))
        if (entry.path().filename().string().find(".tmp.") != std::string::npos)
            stray = true;
    check(!stray, "failure paths remove temporary output");
    check(equilibrium::write_spdata(path, valid).has_value(),
          "a later write succeeds and replaces the target");
}
}  // namespace

int main() {
    try {
        // Formatting must not follow the caller's decimal separator.
        const auto original = std::locale();
        std::locale::global(std::locale(original, new CommaDecimal));
        for (std::size_t lst : {1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 999u})
            check_file(lst);
        std::locale::global(original);
        failures();
    } catch (const std::exception& error) { check(false, error.what()); }
    return equilibrium::test::summary();
}
