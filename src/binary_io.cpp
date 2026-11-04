#include "equilibrium/binary.hpp"
#include "equilibrium/checkpoint.hpp"
#include "equilibrium/validation.hpp"
#include "record_io.hpp"

#include <array>
#include <fstream>
#include <string_view>

namespace equilibrium {
namespace io_detail {

Result<Equilibrium> read_record(Reader& reader,
                                const std::string& path,
                                OutputFormat format,
                                ReadMode mode) {
    RunReport report;
    auto state = reader.read(path, mode == ReadMode::COMPLETE
                                       ? std::optional(std::ref(report))
                                       : std::nullopt);
    if (!state) return Result<Equilibrium>(state.error());
    auto owner = take_equilibrium(std::move(state.value()));
    owner.source.format = format;
    if (mode == ReadMode::COMPLETE) {
        owner.grid = grid_metadata(report.input_params);
        owner.input = std::move(report.input_params);
        owner.provenance = std::move(static_cast<RunProvenance&>(report));
    }
    auto status = validate(as_view(owner));
    if (!status) return Result<Equilibrium>(status.error());
    return owner;
}

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

Result<Equilibrium> read_binary(const std::string& path, ReadMode mode) {
    std::ifstream file(path, std::ios::binary);
    std::array<char, 8> magic{};
    file.read(magic.data(), magic.size());
    if (!file)
        return Result<Equilibrium>("missing or truncated input: " + path);
    file.close();
    if (std::string_view(magic.data(), magic.size()) == "CUMECKP1") {
        InputParams input;
        auto state = read_checkpoint(path, mode == ReadMode::COMPLETE
                                               ? std::optional(std::ref(input))
                                               : std::nullopt);
        if (!state) return Result<Equilibrium>(state.error());
        auto owner = take_equilibrium(std::move(state.value()));
        owner.source.format = OutputFormat::BINARY;
        if (mode == ReadMode::COMPLETE) {
            owner.grid = grid_metadata(input);
            owner.input = std::move(input);
        }
        auto status = validate(as_view(owner));
        if (!status) return Result<Equilibrium>(status.error());
        return owner;
    }
    return io_detail::read_record(*make_binary_reader(), path,
                                  OutputFormat::BINARY, mode);
}

Status write_binary(const std::string& path,
                    const EquilibriumView& view,
                    WriteOptions options) {
    return io_detail::write_record(*make_binary_writer(),
                                   {OutputFormat::BINARY, path}, view, options);
}

}  // namespace equilibrium
