#include "boozer_output_internal.hpp"
#include "equilibrium/boozer_binary.hpp"
#include "equilibrium/boozer_output.hpp"
#include "equilibrium/detail/list_profile.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string>

namespace equilibrium {
namespace {

std::string lower_extension(const std::filesystem::path& path) {
    std::string extension = path.extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char value) {
                       return static_cast<char>(std::tolower(value));
                   });
    return extension;
}

}  // namespace

BoozerOutputSpec resolve_boozer_output_spec(const std::filesystem::path& path) {
    const std::string extension = lower_extension(path);
    if (extension == ".bin") { return {BoozerOutputFormat::BINARY, path}; }
    if (extension == ".nc") { return {BoozerOutputFormat::NETCDF, path}; }
    if (extension == ".h5" || extension == ".hdf5") {
        return {BoozerOutputFormat::HDF5, path};
    }
    throw std::invalid_argument("unrecognized Boozer output suffix '" +
                                extension +
                                "' (expected .bin, .nc, .h5, or .hdf5)");
}

bool boozer_output_format_available(BoozerOutputFormat format) {
    switch (format) {
        case BoozerOutputFormat::BINARY:
            return true;
        case BoozerOutputFormat::NETCDF:
#ifdef EQUILIBRIUM_HAVE_NETCDF
            return true;
#else
            return false;
#endif
        case BoozerOutputFormat::HDF5:
#ifdef EQUILIBRIUM_HAVE_HDF5
            return true;
#else
            return false;
#endif
    }
    return false;
}

std::string_view boozer_output_suffix(BoozerOutputFormat format) {
    switch (format) {
        case BoozerOutputFormat::BINARY:
            return ".bin";
        case BoozerOutputFormat::NETCDF:
            return ".nc";
        case BoozerOutputFormat::HDF5:
            return ".h5";
    }
    return "";
}

void write_boozer_output(const BoozerOutputSpec& spec,
                         const BoozerResult& result,
                         const std::filesystem::path& source_path) {
    if (!boozer_output_format_available(spec.format)) {
        throw std::runtime_error(
            "Boozer output format '" +
            std::string(boozer_output_suffix(spec.format)) +
            "' is not available in this build");
    }
    switch (spec.format) {
        case BoozerOutputFormat::BINARY:
            write_boozer_binary(spec.path, result, source_path);
            return;
        case BoozerOutputFormat::NETCDF:
#ifdef EQUILIBRIUM_HAVE_NETCDF
            boozer_detail::write_boozer_netcdf(spec.path, result, source_path);
            return;
#else
            break;
#endif
        case BoozerOutputFormat::HDF5:
#ifdef EQUILIBRIUM_HAVE_HDF5
            boozer_detail::write_boozer_hdf5(spec.path, result, source_path);
            return;
#else
            break;
#endif
    }
    throw std::runtime_error("unavailable Boozer output format");
}

}  // namespace equilibrium
