#include "equilibrium/binary.hpp"
#include "equilibrium/checkpoint.hpp"
#include "equilibrium/validation.hpp"
#include "read_selection.hpp"
#include "record_io.hpp"

#include <array>
#include <fstream>
#include <string_view>

namespace equilibrium {
namespace io_detail {

Status write_record(Writer& writer,
                    const OutputSpec& output,
                    const EquilibriumView& view,
                    WriteOptions options) {
    auto status = validate(view);
    if (!status) return status;
    if (!view.input || !view.provenance)
        return Status("result output requires saved input and run provenance");
    if (view.profiles && options.profiles != ProfileWritePolicy::OMIT)
        return Status(
            "this schema cannot store solved profiles; explicitly "
            "select ProfileWritePolicy::OMIT to discard them");
    RunReport report;
    static_cast<RunProvenance&>(report) = view.provenance->get();
    report.input_params = view.input->get();
    return writer.write_atomic(as_snapshot_view(view), report, output);
}

}  // namespace io_detail

Result<Equilibrium> read_binary(const std::string& path,
                                const ReadRequest& request) {
    return io_detail::read_native_binary(path, request);
}

Result<BinaryKind> detect_binary_kind(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    std::array<char, 8> bytes{};
    file.read(bytes.data(), bytes.size());
    if (!file) return Result<BinaryKind>("missing or truncated input: " + path);
    const std::string_view magic(bytes.data(), bytes.size());
    if (magic == "CUMES001" || magic == "CUMECKP1") return BinaryKind::NATIVE;
    if (magic.starts_with("MCBOOZ")) return BinaryKind::BOOZER;
    return Result<BinaryKind>("unsupported equilibrium binary: " + path);
}

Status write_binary(const std::string& path,
                    const EquilibriumView& view,
                    WriteOptions options) {
    return io_detail::write_record(*make_binary_writer(),
                                   {OutputFormat::BINARY, path}, view, options);
}

}  // namespace equilibrium
