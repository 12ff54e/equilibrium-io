#include "equilibrium/boozer_binary.hpp"

#include "boozer_output_internal.hpp"
#include "equilibrium/writer_helpers.hpp"

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace equilibrium {
namespace {

constexpr std::array<char, 8> MAGIC{'M', 'C', 'B', 'O', 'O', 'Z', '0', '3'};
constexpr std::int32_t VERSION = 3;
constexpr std::array<char, 8> ASYMMETRIC_MAGIC{'M', 'C', 'B', 'O',
                                               'O', 'Z', '0', '4'};
constexpr std::int32_t ASYMMETRIC_VERSION = 4;
constexpr std::array<char, 8> TRACING_MAGIC{'M', 'C', 'B', 'O',
                                            'O', 'Z', '0', '5'};
constexpr std::int32_t TRACING_VERSION = 5;
constexpr std::array<char, 8> MAGNETIC_MAGIC{'M', 'C', 'B', 'O',
                                             'O', 'Z', '0', '6'};
constexpr std::int32_t MAGNETIC_VERSION = 6;
constexpr std::array<char, 8> PROFILE_MAGIC{'M', 'C', 'B', 'O',
                                            'O', 'Z', '0', '7'};
constexpr std::int32_t PROFILE_VERSION = 7;
constexpr std::array<char, 8> FLUX_MAGIC{'M', 'C', 'B', 'O',
                                         'O', 'Z', '0', '8'};
constexpr std::int32_t FLUX_VERSION = 8;

std::size_t checked_product(std::size_t left, std::size_t right) {
    if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left) {
        throw std::runtime_error("Boozer binary dimension overflow");
    }
    return left * right;
}

class Writer {
   public:
    explicit Writer(const std::filesystem::path& path)
        : path_(path.string()),
          temporary_(io_detail::temp_path_for(path_)),
          output_(std::fopen(temporary_.c_str(), "wbx"), std::fclose) {
        if (!output_) throw std::runtime_error("cannot open Boozer output");
    }
    ~Writer() {
        output_.reset();
        std::remove(temporary_.c_str());
    }
    void bytes(const void* data, std::size_t count) {
        if (std::fwrite(data, 1, count, output_.get()) != count)
            throw std::runtime_error("could not write Boozer output");
    }
    void publish() {
        const auto error =
            io_detail::publish_atomic(output_.release(), temporary_, path_);
        if (!error.empty()) throw std::runtime_error(error);
    }
    template <typename UInt>
    void uint(UInt value) {
        static_assert(std::is_unsigned_v<UInt>);
        std::array<std::uint8_t, sizeof(UInt)> data{};
        for (std::size_t i = 0; i < data.size(); ++i) {
            data[i] = static_cast<std::uint8_t>(value >> (8U * i));
        }
        bytes(data.data(), data.size());
    }
    void i32(std::int32_t value) { uint(std::bit_cast<std::uint32_t>(value)); }
    void f64(double value) { uint(std::bit_cast<std::uint64_t>(value)); }
    void string(std::string_view value) {
        if (value.size() > static_cast<std::size_t>(
                               std::numeric_limits<std::int32_t>::max())) {
            throw std::runtime_error("Boozer output string is too large");
        }
        i32(static_cast<std::int32_t>(value.size()));
        bytes(value.data(), value.size());
    }
    void doubles(std::span<const double> values) {
        for (double value : values) f64(value);
    }
    void integers(std::span<const int> values) {
        for (int value : values) i32(value);
    }

   private:
    std::string path_, temporary_;
    std::unique_ptr<std::FILE, decltype(&std::fclose)> output_;
};

class Reader {
   public:
    explicit Reader(const std::filesystem::path& path)
        : input_(path, std::ios::binary) {
        if (!input_) throw std::runtime_error("cannot open Boozer input");
        remaining_ = std::filesystem::file_size(path);
    }
    void bytes(void* data, std::size_t count) {
        require_extent(count, 1);
        input_.read(static_cast<char*>(data),
                    static_cast<std::streamsize>(count));
        if (!input_) throw std::runtime_error("truncated Boozer input");
        remaining_ -= count;
    }
    template <typename UInt>
    UInt uint() {
        static_assert(std::is_unsigned_v<UInt>);
        std::array<std::uint8_t, sizeof(UInt)> data{};
        bytes(data.data(), data.size());
        UInt value = 0;
        for (std::size_t i = 0; i < data.size(); ++i) {
            value |= static_cast<UInt>(data[i]) << (8U * i);
        }
        return value;
    }
    std::int32_t i32() {
        return std::bit_cast<std::int32_t>(uint<std::uint32_t>());
    }
    double f64() { return std::bit_cast<double>(uint<std::uint64_t>()); }
    std::string string() {
        const auto count = i32();
        if (count < 0 || count > (1 << 24)) {
            throw std::runtime_error("invalid Boozer string length");
        }
        std::string value(static_cast<std::size_t>(count), '\0');
        bytes(value.data(), value.size());
        return value;
    }
    std::vector<double> doubles(std::size_t count) {
        require_extent(count, sizeof(double));
        std::vector<double> values(count);
        for (double& value : values) value = f64();
        return values;
    }
    std::vector<int> integers(std::size_t count) {
        require_extent(count, sizeof(std::int32_t));
        std::vector<int> values(count);
        for (int& value : values) value = i32();
        return values;
    }

   private:
    void require_extent(std::size_t count, std::size_t width) const {
        if (count > remaining_ / width)
            throw std::runtime_error("truncated Boozer input");
    }
    std::ifstream input_;
    std::uintmax_t remaining_ = 0;
};

}  // namespace

void write_boozer_binary(const std::filesystem::path& path,
                         const BoozerResult& result,
                         const std::filesystem::path& source_path) {
    const boozer_detail::RealBoozerSpectrum spectrum =
        boozer_detail::make_real_spectrum(result);
    const bool has_flux_type =
        result.tracing && result.tracing->pflux_type == "list";
    const bool has_magnetic = result.tracing && result.tracing->magnetic;
    Writer writer(path);
    const auto& magic = has_flux_type     ? FLUX_MAGIC
                        : result.profiles ? PROFILE_MAGIC
                        : has_magnetic    ? MAGNETIC_MAGIC
                        : result.tracing  ? TRACING_MAGIC
                        : result.lasym    ? ASYMMETRIC_MAGIC
                                          : MAGIC;
    writer.bytes(magic.data(), magic.size());
    writer.i32(has_flux_type     ? FLUX_VERSION
               : result.profiles ? PROFILE_VERSION
               : has_magnetic    ? MAGNETIC_VERSION
               : result.tracing  ? TRACING_VERSION
               : result.lasym    ? ASYMMETRIC_VERSION
                                 : VERSION);
    writer.string(boozer_detail::COORDINATE_CONVENTION);
    writer.string(result.lasym ? boozer_detail::ASYMMETRIC_FOURIER_CONVENTION
                               : boozer_detail::FOURIER_CONVENTION);
    writer.string(source_path.string());
    writer.i32(result.source_format_version);
    writer.i32(result.source_ns);
    writer.i32(result.source_ntheta);
    writer.i32(result.source_nzeta);
    writer.i32(result.source_mpol);
    writer.i32(result.source_ntor);
    writer.i32(result.nfp);
    writer.i32(result.grid.first_surface);
    writer.i32(result.grid.ntheta);
    writer.i32(result.grid.nzeta);
    writer.i32(result.spectrum.mmax);
    writer.i32(result.spectrum.nmax);
    writer.i32(static_cast<std::int32_t>(result.radial_order));
    if (result.tracing) writer.i32(result.lasym ? 1 : 0);
    writer.f64(result.resonance_tolerance);
    writer.doubles(result.s);
    writer.doubles(result.iota);
    writer.integers(spectrum.m);
    writer.integers(spectrum.n);
    writer.doubles(spectrum.rmncc);
    writer.doubles(spectrum.rmnss);
    writer.doubles(spectrum.zmnsc);
    writer.doubles(spectrum.zmncs);
    writer.doubles(spectrum.numnsc);
    writer.doubles(spectrum.numncs);
    if (result.lasym) {
        writer.doubles(spectrum.rmnsc);
        writer.doubles(spectrum.rmncs);
        writer.doubles(spectrum.zmncc);
        writer.doubles(spectrum.zmnss);
        writer.doubles(spectrum.numncc);
        writer.doubles(spectrum.numnss);
    }
    writer.doubles(result.grid.b);
    writer.doubles(result.grid.sqrtg_b);
    writer.doubles(result.grid.b2j00);
    if (result.tracing) {
        const auto& tracing = *result.tracing;
        writer.string(boozer_detail::TRACING_CONVENTION);
        writer.i32(static_cast<std::int32_t>(tracing.aphi.size()));
        writer.f64(tracing.phiedge);
        writer.doubles(tracing.aphi);
        writer.doubles(tracing.toroidal_flux_derivative);
        writer.doubles(tracing.poloidal_flux_derivative);
        writer.doubles(tracing.axis_r_cos);
        writer.doubles(tracing.axis_r_sin);
        writer.doubles(tracing.axis_z_cos);
        writer.doubles(tracing.axis_z_sin);
        writer.doubles(tracing.axis_half_s);
        for (const auto& field : tracing.axis_half_fields)
            writer.doubles(field);
        if (tracing.magnetic) {
            const auto& magnetic = *tracing.magnetic;
            writer.string(boozer_detail::MAGNETIC_CONVENTION);
            for (const auto& field : magnetic.fields) writer.doubles(field);
            writer.doubles(result.grid.source_theta_shift);
            writer.doubles(result.grid.nu);
        }
    }
    if (result.profiles) {
        writer.string(boozer_detail::PROFILE_CONVENTION);
        writer.doubles(result.profiles->pressure);
        writer.doubles(result.profiles->current_i);
        writer.doubles(result.profiles->current_g);
    }
    if (has_flux_type) writer.string(result.tracing->pflux_type);
    writer.publish();
}

BoozerFile read_boozer_binary(const std::filesystem::path& path) {
    Reader reader(path);
    std::array<char, 8> magic{};
    reader.bytes(magic.data(), magic.size());
    const int version = reader.i32();
    bool lasym = magic == ASYMMETRIC_MAGIC && version == ASYMMETRIC_VERSION;
    const bool has_flux_type = magic == FLUX_MAGIC && version == FLUX_VERSION;
    const bool has_profiles =
        has_flux_type || (magic == PROFILE_MAGIC && version == PROFILE_VERSION);
    const bool has_magnetic = has_profiles || (magic == MAGNETIC_MAGIC &&
                                               version == MAGNETIC_VERSION);
    const bool has_tracing =
        has_magnetic || (magic == TRACING_MAGIC && version == TRACING_VERSION);
    if (!has_tracing && !lasym && (magic != MAGIC || version != VERSION)) {
        throw std::runtime_error("unsupported Boozer binary format");
    }
    if (reader.string() != boozer_detail::COORDINATE_CONVENTION) {
        throw std::runtime_error("unsupported Boozer coordinate convention");
    }
    const auto fourier_convention = reader.string();
    BoozerFile file;
    file.source_path = reader.string();
    auto& result = file.result;
    result.source_format_version = reader.i32();
    result.source_ns = reader.i32();
    result.source_ntheta = reader.i32();
    result.source_nzeta = reader.i32();
    result.source_mpol = reader.i32();
    result.source_ntor = reader.i32();
    result.nfp = reader.i32();
    result.grid.first_surface = reader.i32();
    result.grid.ntheta = reader.i32();
    result.grid.nzeta = reader.i32();
    result.spectrum.mmax = reader.i32();
    result.spectrum.nmax = reader.i32();
    result.radial_order = static_cast<RadialInterpolationOrder>(reader.i32());
    if (has_tracing) {
        const int flag = reader.i32();
        if (flag != 0 && flag != 1)
            throw std::runtime_error("invalid Boozer symmetry flag");
        lasym = flag == 1;
    }
    result.lasym = lasym;
    if (fourier_convention !=
        (lasym ? boozer_detail::ASYMMETRIC_FOURIER_CONVENTION
               : boozer_detail::FOURIER_CONVENTION))
        throw std::runtime_error("unsupported Boozer Fourier convention");
    result.resonance_tolerance = reader.f64();
    result.grid.source_ns = result.source_ns;
    if (result.source_ns < 2 || result.source_ns > (1 << 20) ||
        result.grid.first_surface < 1 ||
        result.grid.first_surface >= result.source_ns)
        throw std::runtime_error("invalid Boozer binary surfaces");
    const int surfaces = result.source_ns - result.grid.first_surface;
    if (result.spectrum.mmax < 0 || result.spectrum.nmax < 0 ||
        result.spectrum.mmax > (1 << 15) || result.spectrum.nmax > (1 << 15)) {
        throw std::runtime_error("invalid Boozer binary mode bounds");
    }
    const std::size_t mode_count =
        checked_product(static_cast<std::size_t>(result.spectrum.mmax) + 1,
                        static_cast<std::size_t>(result.spectrum.nmax) + 1);
    if (surfaces < 1 || result.grid.ntheta < 1 || result.grid.nzeta < 1 ||
        mode_count < 1 || surfaces > (1 << 20) || mode_count > (1 << 20)) {
        throw std::runtime_error("invalid Boozer binary dimensions");
    }
    const std::size_t real_count = checked_product(
        static_cast<std::size_t>(surfaces),
        checked_product(static_cast<std::size_t>(result.grid.ntheta),
                        static_cast<std::size_t>(result.grid.nzeta)));
    const std::size_t coefficient_count =
        checked_product(static_cast<std::size_t>(surfaces), mode_count);
    result.s = reader.doubles(static_cast<std::size_t>(surfaces));
    result.iota = reader.doubles(static_cast<std::size_t>(surfaces));
    boozer_detail::RealBoozerSpectrum spectrum;
    spectrum.lasym = lasym;
    spectrum.mmax = result.spectrum.mmax;
    spectrum.nmax = result.spectrum.nmax;
    spectrum.m = reader.integers(mode_count);
    spectrum.n = reader.integers(mode_count);
    spectrum.rmncc = reader.doubles(coefficient_count);
    spectrum.rmnss = reader.doubles(coefficient_count);
    spectrum.zmnsc = reader.doubles(coefficient_count);
    spectrum.zmncs = reader.doubles(coefficient_count);
    spectrum.numnsc = reader.doubles(coefficient_count);
    spectrum.numncs = reader.doubles(coefficient_count);
    if (lasym) {
        spectrum.rmnsc = reader.doubles(coefficient_count);
        spectrum.rmncs = reader.doubles(coefficient_count);
        spectrum.zmncc = reader.doubles(coefficient_count);
        spectrum.zmnss = reader.doubles(coefficient_count);
        spectrum.numncc = reader.doubles(coefficient_count);
        spectrum.numnss = reader.doubles(coefficient_count);
    }
    result.spectrum = boozer_detail::restore_complex_spectrum(
        spectrum, result.source_ns, result.grid.first_surface);
    result.grid.b = reader.doubles(real_count);
    result.grid.sqrtg_b = reader.doubles(real_count);
    result.grid.b2j00 = reader.doubles(static_cast<std::size_t>(surfaces));
    if (has_tracing) {
        if (reader.string() != boozer_detail::TRACING_CONVENTION)
            throw std::runtime_error("unsupported Boozer tracing convention");
        const int aphi_count = reader.i32();
        const int width = static_cast<int>(result.radial_order);
        if (aphi_count < 1 || aphi_count > (1 << 20) ||
            (width != 2 && width != 4) || result.source_ns <= width ||
            result.source_ntor < 0 || result.source_ntor > (1 << 15) ||
            result.source_ntheta < 1 || result.source_nzeta < 1)
            throw std::runtime_error("invalid Boozer tracing dimensions");
        auto& tracing = result.tracing.emplace();
        tracing.phiedge = reader.f64();
        tracing.aphi = reader.doubles(static_cast<std::size_t>(aphi_count));
        tracing.toroidal_flux_derivative =
            reader.doubles(static_cast<std::size_t>(result.source_ns));
        tracing.poloidal_flux_derivative =
            reader.doubles(static_cast<std::size_t>(result.source_ns));
        const auto axis_modes =
            static_cast<std::size_t>(result.source_ntor) + 1;
        tracing.axis_r_cos = reader.doubles(axis_modes);
        tracing.axis_r_sin = reader.doubles(axis_modes);
        tracing.axis_z_cos = reader.doubles(axis_modes);
        tracing.axis_z_sin = reader.doubles(axis_modes);
        tracing.axis_half_s = reader.doubles(static_cast<std::size_t>(width));
        const auto half_count = checked_product(
            static_cast<std::size_t>(width),
            checked_product(static_cast<std::size_t>(result.source_ntheta),
                            static_cast<std::size_t>(result.source_nzeta)));
        for (auto& field : tracing.axis_half_fields)
            field = reader.doubles(half_count);
        if (has_magnetic) {
            if (reader.string() != boozer_detail::MAGNETIC_CONVENTION)
                throw std::runtime_error(
                    "unsupported Boozer magnetic convention");
            auto& magnetic = tracing.magnetic.emplace();
            for (auto& field : magnetic.fields)
                field = reader.doubles(real_count);
            result.grid.source_theta_shift = reader.doubles(real_count);
            result.grid.nu = reader.doubles(real_count);
        }
    }
    if (has_profiles) {
        if (reader.string() != boozer_detail::PROFILE_CONVENTION)
            throw std::runtime_error(
                "unsupported Boozer surface profile convention");
        auto& profiles = result.profiles.emplace();
        profiles.pressure = reader.doubles(static_cast<std::size_t>(surfaces));
        profiles.current_i = reader.doubles(static_cast<std::size_t>(surfaces));
        profiles.current_g = reader.doubles(static_cast<std::size_t>(surfaces));
    }
    if (has_flux_type) result.tracing->pflux_type = reader.string();
    validate_boozer_result(result);
    return file;
}

}  // namespace equilibrium
