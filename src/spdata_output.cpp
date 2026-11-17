#include "equilibrium/spdata_output.hpp"

#include "equilibrium/spdata_validation.hpp"
#include "equilibrium/writer_helpers.hpp"

#include <array>
#include <cstdio>
#include <iomanip>
#include <locale>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace equilibrium {
namespace {

struct CloseFile {
    void operator()(std::FILE* file) const { std::fclose(file); }
};

class SpdataWriter {
   public:
    explicit SpdataWriter(const std::string& path)
        : path_(path),
          temporary_(io_detail::temp_path_for(path)),
          file_(std::fopen(temporary_.c_str(), "wbx")) {
        if (!file_) throw std::runtime_error("cannot open spdata output");
    }
    ~SpdataWriter() {
        file_.reset();
        std::remove(temporary_.c_str());
    }
    void write(std::string_view bytes) {
        if (std::fwrite(bytes.data(), 1, bytes.size(), file_.get()) !=
            bytes.size())
            throw std::runtime_error("cannot write spdata output");
    }
    Status publish() {
        const auto error =
            io_detail::publish_atomic(file_.release(), temporary_, path_);
        if (!error.empty()) return error;
        return {};
    }

   private:
    std::string path_, temporary_;
    std::unique_ptr<std::FILE, CloseFile> file_;
};

}  // namespace

Status write_spdata(const std::string& path, const SpdataView& data) {
    if (auto status = validate_spdata(data); !status) return status;
    if (path.empty() || path.find('\0') != std::string::npos)
        return std::string("invalid spdata output path");
    try {
        SpdataWriter writer(path);
        std::ostringstream output;
        output.imbue(std::locale::classic());
        output.exceptions(std::ios::badbit | std::ios::failbit);
        output << data.meta_data << '\n';
        output << std::setw(4) << data.lsp << std::setw(4) << data.lst
               << std::setw(4) << 4 << std::setw(4) << 8 << '\n';
        output << std::scientific << std::uppercase << std::setprecision(10);
        output << std::setw(18) << data.psi_wall << std::setw(18)
               << data.psi_sep << '\n';
        const auto flush = [&] {
            writer.write(output.view());
            output.str(std::string{});
        };
        flush();

        const std::array fields{data.b, data.r, data.z, data.j};
        const std::array profiles{data.q, data.f,       data.i,
                                  data.p, data.r_minor, data.psi_t};
        // Buffer one radial record, not the whole output file.
        for (std::size_t radial = 0; radial < data.lsp; ++radial) {
            for (const auto field : fields) {
                for (std::size_t c = 0; c < 9; ++c) {
                    for (std::size_t t = 0; t <= data.lst; ++t) {
                        output
                            << std::setw(18)
                            << field[data.field_index(radial, c, t % data.lst)];
                        if (t % 4 == 3) output << '\n';
                    }
                    if ((data.lst + 1) % 4 != 0) output << '\n';
                }
            }
            for (const auto profile : profiles) {
                for (std::size_t c = 0; c < 3; ++c)
                    output << std::setw(18)
                           << profile[data.profile_index(radial, c)];
                output << '\n';
            }
            flush();
        }
        // Axisymmetric output has no ripple. Use the supplied axis radius,
        // which need not be one when the coefficients carry physical units.
        output << std::setw(4) << 0 << std::setw(4) << 0 << '\n';
        output << std::setw(18) << data.r.front() << std::setw(18) << 0.0
               << std::setw(18) << 0.0 << '\n';
        output << std::setw(18) << 0.0 << std::setw(18) << 0.0 << '\n';
        flush();
        return writer.publish();
    } catch (const std::exception& error) {
        return std::string("spdata output: ") + error.what();
    }
}

}  // namespace equilibrium
