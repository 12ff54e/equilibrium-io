#include "equilibrium/miller_batch.hpp"

#include <charconv>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {
template <class T>
T number(std::string_view text) {
    T value{};
    const auto [end, error] =
        std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size())
        throw std::invalid_argument("invalid option value: " +
                                    std::string(text));
    return value;
}
equilibrium::MillerRange parse_range(std::string_view text) {
    const auto first = text.find(':');
    if (first == std::string_view::npos) {
        const double fixed = number<double>(text);
        return {fixed, fixed, 1.0};
    }
    const auto second = text.find(':', first + 1);
    if (second == std::string_view::npos)
        throw std::invalid_argument("expected VALUE or FIRST:LAST:STEP");
    return {number<double>(text.substr(0, first)),
            number<double>(text.substr(first + 1, second - first - 1)),
            number<double>(text.substr(second + 1))};
}
void print_shape(const equilibrium::LcfsShape& shape, int cocos) {
    std::cout << std::setprecision(17) << "{\n  \"cocos\":" << cocos
              << ",\n  \"r_center\":" << shape.r_center
              << ",\n  \"minor_radius\":" << shape.minor_radius
              << ",\n  \"z_reference\":" << shape.z_reference
              << ",\n  \"elongation_upper\":" << shape.elongation_upper
              << ",\n  \"elongation_lower\":" << shape.elongation_lower
              << ",\n  \"triangularity_upper\":" << shape.triangularity_upper
              << ",\n  \"triangularity_lower\":" << shape.triangularity_lower
              << "\n}\n";
}
}  // namespace
int main(int argc, char** argv) {
    try {
        equilibrium::GfileInputSettings settings;
        std::string input, output, batch_directory;
        equilibrium::MillerScan scan;
        bool batch_option = false;
        bool shape_only = false;
        for (int arg = 1; arg < argc; ++arg) {
            const std::string_view token(argv[arg]);
            const auto next = [&]() -> std::string_view {
                if (++arg == argc)
                    throw std::invalid_argument(std::string(token) +
                                                " requires a value");
                return argv[arg];
            };
            if (token == "--output" || token == "-o")
                output = next();
            else if (token == "--cocos")
                settings.cocos = number<int>(next());
            else if (token == "--mpol")
                settings.mpol = number<int>(next());
            else if (token == "--ntheta")
                settings.ntheta = number<int>(next());
            else if (token == "--psi-ratio")
                settings.psi_ratio = number<double>(next());
            else if (token == "--profile-points")
                settings.profile_points = number<std::size_t>(next());
            else if (token == "--boundary-points")
                settings.boundary_points = number<std::size_t>(next());
            else if (token == "--boundary-tolerance")
                settings.boundary_tolerance = number<double>(next());
            else if (token == "--miller-batch")
                batch_directory = next();
            else if (token == "--max-cases") {
                batch_option = true;
                scan.max_cases = number<std::size_t>(next());
            } else if (token == "--r0" || token == "--kappa-upper" ||
                       token == "--kappa-lower" || token == "--delta-upper" ||
                       token == "--delta-lower") {
                batch_option = true;
                const auto range = parse_range(next());
                if (token == "--r0")
                    scan.r0 = range;
                else if (token == "--kappa-upper")
                    scan.kappa_upper = range;
                else if (token == "--kappa-lower")
                    scan.kappa_lower = range;
                else if (token == "--delta-upper")
                    scan.delta_upper = range;
                else
                    scan.delta_lower = range;
            } else if (token == "--shape")
                shape_only = true;
            else if (token == "--help" || token == "-h") {
                std::cout
                    << "usage: equilibrium-gfile INPUT --output input.json "
                       "[--cocos N] [--mpol N]\n"
                       "       [--ntheta N] [--psi-ratio X] [--profile-points "
                       "N] [--boundary-points N]\n"
                       "       equilibrium-gfile INPUT --shape [--cocos N]\n"
                       "       equilibrium-gfile INPUT --miller-batch NEW_DIR "
                       "[--r0 RANGE]\n"
                       "          [--kappa-upper RANGE] [--kappa-lower RANGE]\n"
                       "          [--delta-upper RANGE] [--delta-lower RANGE]\n"
                       "          [--mpol N] [--boundary-tolerance X] "
                       "[--max-cases N]\n"
                       "       RANGE is VALUE or FIRST:LAST:STEP. Omitted "
                       "ranges use source LCFS values.\n";
                return 0;
            } else if (!token.empty() && token.front() == '-')
                throw std::invalid_argument("unknown option: " +
                                            std::string(token));
            else if (input.empty())
                input = token;
            else
                throw std::invalid_argument(
                    "only one G-file input is accepted");
        }
        if (input.empty())
            throw std::invalid_argument("a G-file input is required");
        if (!batch_directory.empty()) {
            if (shape_only || !output.empty())
                throw std::invalid_argument(
                    "--miller-batch requires its own output directory");
            const auto count = equilibrium::write_miller_batch(
                input, batch_directory, scan, settings);
            std::cout << "wrote " << count << " Miller inputs in "
                      << batch_directory << '\n';
            return 0;
        }
        if (batch_option)
            throw std::invalid_argument(
                "scan ranges and --max-cases require --miller-batch");
        if (shape_only) {
            if (!output.empty())
                throw std::invalid_argument(
                    "--shape cannot be combined with --output");
            const auto lcfs =
                equilibrium::read_gfile_lcfs(input, settings.cocos);
            print_shape(equilibrium::calculate_lcfs_shape(lcfs), lcfs.cocos);
            return 0;
        }
        if (output.empty())
            throw std::invalid_argument(
                "--output is required for input conversion");
        const auto converted =
            equilibrium::gfile_to_cumes_input(input, settings);
        const auto status = equilibrium::write_gfile_input(output, converted);
        if (!status) throw std::runtime_error(status.error());
        std::cout << "wrote " << output
                  << "; source COCOS=" << converted.lcfs.cocos
                  << "; boundary fit error=" << converted.boundary_max_error
                  << " m"
                  << "; relative poloidal-flux error="
                  << converted.poloidal_flux_relative_error << '\n';
    } catch (const std::exception& error) {
        std::cerr << "equilibrium-gfile: " << error.what() << '\n';
        return 1;
    }
}
