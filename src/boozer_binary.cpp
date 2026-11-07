#include "equilibrium/boozer_binary.hpp"

#include "boozer_output_internal.hpp"
#include "equilibrium/boozer_validation.hpp"
#include "equilibrium/writer_helpers.hpp"
#include "read_selection.hpp"

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
#include <utility>
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

class FramingError : public std::runtime_error {
   public:
    using std::runtime_error::runtime_error;
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
        if (!input_) throw FramingError("truncated Boozer input");
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
            throw FramingError("invalid Boozer string length");
        }
        std::string value(static_cast<std::size_t>(count), '\0');
        bytes(value.data(), value.size());
        return value;
    }
    void skip(std::size_t count, std::size_t width = 1) {
        require_extent(count, width);
        const auto bytes = count * width;
        if (!std::in_range<std::streamoff>(bytes))
            throw FramingError("Boozer extent exceeds stream range");
        input_.seekg(static_cast<std::streamoff>(bytes), std::ios::cur);
        if (!input_) throw FramingError("cannot skip Boozer payload");
        remaining_ -= bytes;
    }
    std::vector<double> doubles(std::size_t count, bool load = true) {
        if (!load) {
            skip(count, sizeof(double));
            return {};
        }
        require_extent(count, sizeof(double));
        std::vector<double> values(count);
        for (double& value : values) value = f64();
        return values;
    }
    std::vector<int> integers(std::size_t count, bool load = true) {
        if (!load) {
            skip(count, sizeof(std::int32_t));
            return {};
        }
        require_extent(count, sizeof(std::int32_t));
        std::vector<int> values(count);
        for (int& value : values) value = i32();
        return values;
    }

   private:
    void require_extent(std::size_t count, std::size_t width) const {
        if (count > remaining_ / width)
            throw FramingError("truncated Boozer input");
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
    writer.i32(spectrum.mmax);
    writer.i32(spectrum.nmax);
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

namespace {
BoozerFile read_selected(const std::filesystem::path& path,
                         const BoozerReadOptions& options,
                         bool strict) {
    io_detail::ReadSelection selection(options.groups, strict);
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
    using Q = Quantity;
    if (options.spectrum != SpectrumStorage::REAL &&
        options.spectrum != SpectrumStorage::COMPLEX &&
        options.spectrum != SpectrumStorage::BOTH)
        throw std::invalid_argument("unsupported Boozer spectrum storage");
    result.s = reader.doubles(surfaces);
    const auto coordinates = validate_boozer_coordinates(as_view(result));
    if (!coordinates) throw std::runtime_error(coordinates.error());
    std::array<std::string, static_cast<std::size_t>(Q::COUNT)> errors;
    std::string blocked;
    const auto read_block = [&](Q q, auto decode) {
        if (!blocked.empty()) {
            errors[static_cast<std::size_t>(q)] = blocked;
            return;
        }
        try {
            decode();
        } catch (const FramingError& e) {
            if (strict) throw;
            blocked = errors[static_cast<std::size_t>(q)] = e.what();
        }
    };
    const auto doubles = [&](Q q, auto& values, std::size_t count) {
        read_block(q,
                   [&] { values = reader.doubles(count, selection.wants(q)); });
    };
    const auto real_families = [](auto& r) {
        return std::array{&r.rmncc,  &r.rmnss,  &r.zmnsc,  &r.zmncs,
                          &r.numnsc, &r.numncs, &r.rmnsc,  &r.rmncs,
                          &r.zmncc,  &r.zmnss,  &r.numncc, &r.numnss};
    };
    const auto drop = [&](Q q) {
        if (q == Q::ROTATIONAL_TRANSFORM) result.iota.clear();
        if (q == Q::GEOMETRY || q == Q::SPECTRAL_SHIFT) {
            if (result.real_spectrum) {
                auto families = real_families(*result.real_spectrum);
                for (std::size_t i = 0; i < families.size(); ++i)
                    if ((i % 6 < 4) == (q == Q::GEOMETRY)) families[i]->clear();
            }
            if (q == Q::GEOMETRY) {
                result.spectrum.r.clear();
                result.spectrum.z.clear();
            } else
                result.spectrum.nu.clear();
        }
        if (q == Q::FIELD_STRENGTH) result.grid.b.clear();
        if (q == Q::STRENGTH_JACOBIAN) result.grid.sqrtg_b.clear();
        if (q == Q::B2_JACOBIAN) result.grid.b2j00.clear();
        if (q == Q::TOROIDAL_SHIFT) result.grid.nu.clear();
        if (q == Q::SOURCE_THETA_SHIFT) result.grid.source_theta_shift.clear();
        if (result.profiles) {
            if (q == Q::PRESSURE) result.profiles->pressure.clear();
            if (q == Q::SURFACE_CURRENTS) {
                result.profiles->current_i.clear();
                result.profiles->current_g.clear();
            }
        }
        if (result.tracing) {
            auto& t = *result.tracing;
            if (q == Q::FLUX_INPUT) {
                t.phiedge = 0;
                t.aphi.clear();
                t.pflux_type.clear();
            }
            if (q == Q::FLUX_DERIVATIVES) {
                t.toroidal_flux_derivative.clear();
                t.poloidal_flux_derivative.clear();
            }
            if (q == Q::AXIS_GEOMETRY) {
                t.axis_r_cos.clear();
                t.axis_r_sin.clear();
                t.axis_z_cos.clear();
                t.axis_z_sin.clear();
            }
            if (q == Q::AXIS_FIELDS) {
                t.axis_half_s.clear();
                for (auto& f : t.axis_half_fields) f.clear();
            }
        }
        if (result.independent_magnetic)
            for (std::size_t i = 0;
                 i < result.independent_magnetic->fields.size(); ++i)
                if (half_field_quantity(i) == q)
                    result.independent_magnetic->fields[i].clear();
    };
    const auto check = [&](Q q, bool present) {
        std::string reason = errors[static_cast<std::size_t>(q)];
        if (selection.wants(q) && present && reason.empty()) {
            const auto status = validate_quantity(as_view(result), q);
            if (!status) reason = status.error();
        }
        if (!selection.accept(q, present, reason)) drop(q);
    };
    selection.accept(Q::GRID, true);
    doubles(Q::ROTATIONAL_TRANSFORM, result.iota, surfaces);
    check(Q::ROTATIONAL_TRANSFORM, true);

    RealBoozerSpectrum spectrum;
    spectrum.lasym = lasym;
    spectrum.mmax = result.spectrum.mmax;
    spectrum.nmax = result.spectrum.nmax;
    const bool spectra = selection.wants(DataGroup::GEOMETRY);
    read_block(Q::GEOMETRY, [&] {
        spectrum.m = reader.integers(mode_count, spectra);
        spectrum.n = reader.integers(mode_count, spectra);
    });
    if (!blocked.empty())
        errors[static_cast<std::size_t>(Q::SPECTRAL_SHIFT)] = blocked;
    auto families = real_families(spectrum);
    for (std::size_t i = 0; i < (lasym ? 12U : 6U); ++i)
        doubles(i % 6 < 4 ? Q::GEOMETRY : Q::SPECTRAL_SHIFT, *families[i],
                coefficient_count);
    if (spectra) {
        result.real_spectrum = std::move(spectrum);
        if (fourier_convention !=
            (lasym ? boozer_detail::ASYMMETRIC_FOURIER_CONVENTION
                   : boozer_detail::FOURIER_CONVENTION))
            for (auto q : group_quantities(DataGroup::GEOMETRY))
                errors[static_cast<std::size_t>(q)] =
                    "unsupported Boozer Fourier convention";
    }
    for (auto q : group_quantities(DataGroup::GEOMETRY)) check(q, true);
    const bool geometry = result.has_quantity(Q::GEOMETRY);
    const bool shift = result.has_quantity(Q::SPECTRAL_SHIFT);
    if (geometry || shift) {
        if (options.spectrum != SpectrumStorage::REAL)
            result.spectrum = boozer_detail::restore_complex_spectrum(
                *result.real_spectrum, result.source_ns,
                result.grid.first_surface, geometry, shift);
        if (options.spectrum == SpectrumStorage::COMPLEX)
            result.real_spectrum.reset();
    } else
        result.real_spectrum.reset();

    doubles(Q::FIELD_STRENGTH, result.grid.b, real_count);
    doubles(Q::STRENGTH_JACOBIAN, result.grid.sqrtg_b, real_count);
    doubles(Q::B2_JACOBIAN, result.grid.b2j00, surfaces);
    for (auto q : group_quantities(DataGroup::STRENGTH)) check(q, true);

    if (has_tracing) {
        BoozerTracingData t;
        std::string convention;
        int aphi_count = 0;
        read_block(Q::FLUX_INPUT, [&] {
            convention = reader.string();
            aphi_count = reader.i32();
            if (aphi_count < 0 || aphi_count > (1 << 20))
                throw FramingError("invalid flux-input extent");
            const double phiedge = reader.f64();
            if (selection.wants(Q::FLUX_INPUT)) t.phiedge = phiedge;
        });
        doubles(Q::FLUX_INPUT, t.aphi, aphi_count);
        doubles(Q::FLUX_DERIVATIVES, t.toroidal_flux_derivative,
                result.source_ns);
        doubles(Q::FLUX_DERIVATIVES, t.poloidal_flux_derivative,
                result.source_ns);
        std::size_t axis_modes = 0;
        read_block(Q::AXIS_GEOMETRY, [&] {
            if (result.source_ntor < 0 || result.source_ntor > (1 << 15))
                throw FramingError("invalid axis mode count");
            axis_modes = static_cast<std::size_t>(result.source_ntor) + 1;
        });
        for (auto* values :
             {&t.axis_r_cos, &t.axis_r_sin, &t.axis_z_cos, &t.axis_z_sin})
            doubles(Q::AXIS_GEOMETRY, *values, axis_modes);
        std::size_t width = 0, half_count = 0;
        read_block(Q::AXIS_FIELDS, [&] {
            const int order = static_cast<int>(result.radial_order);
            if (order < 0 || order > (1 << 20) || result.source_ntheta < 1 ||
                result.source_nzeta < 1)
                throw FramingError("invalid axis field extents");
            width = order;
            half_count = checked_product(
                width,
                checked_product(result.source_ntheta, result.source_nzeta));
        });
        doubles(Q::AXIS_FIELDS, t.axis_half_s, width);
        for (auto& values : t.axis_half_fields)
            doubles(Q::AXIS_FIELDS, values, half_count);
        if (selection.wants(DataGroup::TRACING)) result.tracing = std::move(t);
        if (!convention.empty() &&
            convention != boozer_detail::TRACING_CONVENTION)
            for (auto q : group_quantities(DataGroup::TRACING))
                errors[static_cast<std::size_t>(q)] =
                    "unsupported Boozer tracing convention";
    }
    if (has_magnetic) {
        std::string convention;
        read_block(Q::MAGNETIC_JACOBIAN, [&] { convention = reader.string(); });
        BoozerMagneticData magnetic;
        for (std::size_t i = 0; i < magnetic.fields.size(); ++i)
            doubles(half_field_quantity(i), magnetic.fields[i], real_count);
        doubles(Q::SOURCE_THETA_SHIFT, result.grid.source_theta_shift,
                real_count);
        doubles(Q::TOROIDAL_SHIFT, result.grid.nu, real_count);
        if (selection.wants(Q::MAGNETIC_JACOBIAN) ||
            selection.wants(Q::CONTRAVARIANT_FIELD) ||
            selection.wants(Q::COVARIANT_FIELD))
            result.independent_magnetic = std::move(magnetic);
        if (!convention.empty() &&
            convention != boozer_detail::MAGNETIC_CONVENTION)
            for (auto q : group_quantities(DataGroup::MAGNETIC))
                errors[static_cast<std::size_t>(q)] =
                    "unsupported Boozer magnetic convention";
    }
    if (has_profiles) {
        std::string convention;
        read_block(Q::PRESSURE, [&] { convention = reader.string(); });
        BoozerSurfaceProfiles profiles;
        doubles(Q::PRESSURE, profiles.pressure, surfaces);
        doubles(Q::SURFACE_CURRENTS, profiles.current_i, surfaces);
        doubles(Q::SURFACE_CURRENTS, profiles.current_g, surfaces);
        if (selection.wants(DataGroup::PROFILES))
            result.profiles = std::move(profiles);
        if (!convention.empty() &&
            convention != boozer_detail::PROFILE_CONVENTION)
            for (auto q : group_quantities(DataGroup::PROFILES))
                errors[static_cast<std::size_t>(q)] =
                    "unsupported Boozer profile convention";
    }
    if (has_flux_type && selection.wants(Q::FLUX_INPUT))
        read_block(Q::FLUX_INPUT, [&] {
            const auto type = reader.string();
            if (result.tracing) result.tracing->pflux_type = type;
        });
    for (auto q : group_quantities(DataGroup::TRACING)) check(q, has_tracing);
    for (auto q : group_quantities(DataGroup::MAGNETIC)) check(q, has_magnetic);
    for (auto q : group_quantities(DataGroup::PROFILES)) check(q, has_profiles);

    // Validate relations after local checks, so an invalid optional value does
    // not invalidate an independent required quantity. Resolve required-first.
    for (bool required_pass : {true, false})
        for (auto relation :
             {BoozerConsistency::FLUX_IOTA, BoozerConsistency::SURFACE_CURRENTS,
              BoozerConsistency::MAGNETIC}) {
            const auto quantities = consistency_quantities(relation);
            bool available = true, required = false;
            for (auto q : quantities) {
                available = available && selection.available(q);
                required = required ||
                           selection.requirement(q) == Requirement::REQUIRED;
            }
            if (!available || required != required_pass) continue;
            const auto status = validate_consistency(as_view(result), relation);
            if (!status) {
                if (strict) throw std::invalid_argument(status.error());
                selection.conflict(quantities, status.error(), drop);
            }
        }
    if (result.profiles && !result.has_quantity(Q::PRESSURE) &&
        !result.has_quantity(Q::SURFACE_CURRENTS))
        result.profiles.reset();
    if (result.tracing) {
        bool any = false;
        for (auto q : group_quantities(DataGroup::TRACING))
            any = any || result.has_quantity(q);
        if (!any) result.tracing.reset();
    }
    if (result.independent_magnetic) {
        bool any = false;
        for (const auto& f : result.independent_magnetic->fields)
            any = any || !f.empty();
        if (!any) result.independent_magnetic.reset();
    }
    if (result.tracing && result.independent_magnetic) {
        result.tracing->magnetic = std::move(result.independent_magnetic);
        result.independent_magnetic.reset();
    }
    selection.finish();
    if (strict) validate_boozer_result(result);
    result.read_diagnostics = std::move(selection.diagnostics);
    return file;
}
}  // namespace

BoozerFile read_boozer_binary(const std::filesystem::path& path) {
    return read_selected(path, BoozerReadOptions{}, true);
}

BoozerFile read_boozer_binary(const std::filesystem::path& path,
                              const BoozerReadOptions& options) {
    return read_selected(path, options, false);
}

}  // namespace equilibrium
