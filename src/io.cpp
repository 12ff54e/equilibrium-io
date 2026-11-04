#include "equilibrium/io.hpp"

#include "record_io.hpp"

#include <array>
#include <fstream>
#include <string_view>

namespace equilibrium {

Result<OutputFormat> detect_format(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    std::array<char, 8> bytes{};
    file.read(bytes.data(), bytes.size());
    if (!file)
        return Result<OutputFormat>("missing or truncated input: " + path);
    const std::string_view magic(bytes.data(), bytes.size());
    if (magic == "CUMES001" || magic == "CUMECKP1") return OutputFormat::BINARY;
    if (magic.substr(0, 3) == "CDF" &&
        (bytes[3] == 1 || bytes[3] == 2 || bytes[3] == 5))
        return OutputFormat::NETCDF;
    if (magic == std::string_view("\211HDF\r\n\032\n", 8)) {
        // NetCDF-4 shares HDF5's signature. The container suffix selects its
        // dimension conventions; callers can also select the codec explicitly.
        const auto spec = resolve_output_spec(path);
        if (spec && spec.value().format == OutputFormat::NETCDF)
            return OutputFormat::NETCDF;
        return OutputFormat::HDF5;
    }
    return Result<OutputFormat>("unsupported equilibrium container: " + path);
}

Result<Equilibrium> read(const std::string& path, ReadMode mode) {
    auto format = detect_format(path);
    if (!format) return Result<Equilibrium>(format.error());
    return read(path, format.value(), mode);
}

Result<Equilibrium> read(const std::string& path,
                         OutputFormat format,
                         ReadMode mode) {
    if (format == OutputFormat::BINARY) return read_binary(path, mode);
    auto reader = make_reader(format);
    if (!reader) return Result<Equilibrium>("input codec is unavailable");
    return io_detail::read_record(*reader, path, format, mode);
}

Status write(const std::string& path,
             const EquilibriumView& view,
             WriteOptions options) {
    auto spec = resolve_output_spec(path);
    if (!spec) return Status(spec.error());
    auto writer = create_writer(spec.value().format);
    if (!writer) return Status("output codec is unavailable");
    return io_detail::write_record(*writer, spec.value(), view, options);
}

Status convert(const std::string& input, const std::string& output) {
    auto owner = read(input);
    if (!owner) return Status(owner.error());
    return write(output, owner.value());
}

}  // namespace equilibrium
