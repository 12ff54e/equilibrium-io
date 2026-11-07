#include "boozer_output_internal.hpp"
#include "equilibrium/boozer_binary.hpp"
#include "equilibrium/boozer_output.hpp"
#include "equilibrium/boozer_validation.hpp"
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
std::size_t checked_product(std::size_t left, std::size_t right) {
    if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left) {
        throw std::invalid_argument("Boozer result dimension overflow");
    }
    return left * right;
}

}  // namespace

namespace boozer_detail {

void validate_real_spectrum(const RealBoozerSpectrum& spectrum,
                            int surfaces,
                            bool geometry,
                            bool shift) {
    if (surfaces < 1 || spectrum.mmax < 0 || spectrum.nmax < 0)
        throw std::invalid_argument("invalid real Boozer dimensions");
    const auto modes =
        checked_product(static_cast<std::size_t>(spectrum.mmax) + 1,
                        static_cast<std::size_t>(spectrum.nmax) + 1);
    if (modes > (1U << 20) || spectrum.m.size() != modes ||
        spectrum.n.size() != modes)
        throw std::invalid_argument("invalid real Boozer modes");
    const auto count =
        checked_product(modes, static_cast<std::size_t>(surfaces));
    const std::array families{
        &spectrum.rmncc,  &spectrum.rmnss,  &spectrum.zmnsc,  &spectrum.zmncs,
        &spectrum.numnsc, &spectrum.numncs, &spectrum.rmnsc,  &spectrum.rmncs,
        &spectrum.zmncc,  &spectrum.zmnss,  &spectrum.numncc, &spectrum.numnss};
    for (std::size_t i = 0; i < families.size(); ++i) {
        if (i % 6 < 4 ? !geometry : !shift) continue;
        const auto& family = *families[i];
        if (family.size() != (i < 6 || spectrum.lasym ? count : 0) ||
            !std::ranges::all_of(family,
                                 [](double x) { return std::isfinite(x); }))
            throw std::invalid_argument(
                "invalid real Boozer coefficient array");
    }
    for (std::size_t k = 0; k < modes; ++k)
        if (spectrum.m[k] != static_cast<int>(k % (spectrum.mmax + 1)) ||
            spectrum.n[k] != static_cast<int>(k / (spectrum.mmax + 1)))
            throw std::invalid_argument("noncanonical real Boozer modes");
}

void validate_group(const BoozerResult& result, DataGroup group) {
    const auto view = as_view(result);
    auto status = validate_boozer_coordinates(view);
    if (!status) throw std::invalid_argument(status.error());
    for (auto q : group_quantities(group)) {
        if ((group == DataGroup::TRACING && !result.tracing) ||
            (group == DataGroup::MAGNETIC && !result.magnetic_data()) ||
            (group == DataGroup::PROFILES && !result.profiles))
            continue;
        status = validate_quantity(view, q);
        if (!status) throw std::invalid_argument(status.error());
    }
    if (group == DataGroup::TRACING || group == DataGroup::PROFILES ||
        group == DataGroup::MAGNETIC) {
        const auto check =
            group == DataGroup::TRACING    ? BoozerConsistency::FLUX_IOTA
            : group == DataGroup::PROFILES ? BoozerConsistency::SURFACE_CURRENTS
                                           : BoozerConsistency::MAGNETIC;
        status = validate_consistency(view, check);
        if (!status) throw std::invalid_argument(status.error());
    }
}

void validate_result(const BoozerResult& result) {
    if (result.source_format_version < 1 || !(result.resonance_tolerance > 0) ||
        !std::isfinite(result.resonance_tolerance))
        throw std::invalid_argument("invalid Boozer result metadata");
    const auto iota_status =
        validate_quantity(as_view(result), Quantity::ROTATIONAL_TRANSFORM);
    if (!iota_status) throw std::invalid_argument(iota_status.error());
    for (auto group :
         {DataGroup::GEOMETRY, DataGroup::STRENGTH, DataGroup::TRACING,
          DataGroup::MAGNETIC, DataGroup::PROFILES})
        validate_group(result, group);
    // Current output schemas store these records together. Selected reads may
    // expose them independently, but cannot be written as a complete file.
    if (result.independent_magnetic ||
        (result.profiles && (!result.tracing || !result.tracing->magnetic)) ||
        (result.tracing && result.tracing->pflux_type == "list" &&
         !result.profiles))
        throw std::invalid_argument(
            "selected Boozer groups cannot form a complete output");
}

RealBoozerSpectrum make_real_spectrum(const BoozerResult& result) {
    validate_result(result);
    if (result.real_spectrum) {
        validate_real_spectrum(*result.real_spectrum,
                               result.source_ns - result.grid.first_surface);
        return *result.real_spectrum;
    }
    RealBoozerSpectrum output;
    output.lasym = result.lasym;
    output.mmax = result.spectrum.mmax;
    output.nmax = result.spectrum.nmax;
    const int real_modes = (output.mmax + 1) * (output.nmax + 1);
    const int signed_m_count = 2 * output.mmax + 1;
    const int signed_mode_count = signed_m_count * (2 * output.nmax + 1);
    const int surfaces = result.source_ns - result.grid.first_surface;
    output.m.reserve(static_cast<std::size_t>(real_modes));
    output.n.reserve(static_cast<std::size_t>(real_modes));
    for (int n = 0; n <= output.nmax; ++n) {
        for (int m = 0; m <= output.mmax; ++m) {
            output.m.push_back(m);
            output.n.push_back(n);
        }
    }
    const std::size_t coefficient_count = static_cast<std::size_t>(surfaces) *
                                          static_cast<std::size_t>(real_modes);
    output.rmncc.resize(coefficient_count);
    output.rmnss.resize(coefficient_count);
    output.zmnsc.resize(coefficient_count);
    output.zmncs.resize(coefficient_count);
    output.numnsc.resize(coefficient_count);
    output.numncs.resize(coefficient_count);
    if (output.lasym) {
        output.rmnsc.resize(coefficient_count);
        output.rmncs.resize(coefficient_count);
        output.zmncc.resize(coefficient_count);
        output.zmnss.resize(coefficient_count);
        output.numncc.resize(coefficient_count);
        output.numnss.resize(coefficient_count);
    }

    const auto signed_mode = [&](int m, int n) {
        return (n + output.nmax) * signed_m_count + (m + output.mmax);
    };
    constexpr double FOUR_PI_SQUARED =
        4.0 * std::numbers::pi * std::numbers::pi;
    for (int surface = 0; surface < surfaces; ++surface) {
        for (int mode = 0; mode < real_modes; ++mode) {
            const int m = output.m[static_cast<std::size_t>(mode)];
            const int n = output.n[static_cast<std::size_t>(mode)];
            const std::size_t plus = static_cast<std::size_t>(
                surface * signed_mode_count + signed_mode(m, n));
            const std::size_t minus = static_cast<std::size_t>(
                surface * signed_mode_count + signed_mode(m, -n));
            const std::size_t target =
                static_cast<std::size_t>(surface * real_modes + mode);
            const double mode_weight = static_cast<double>(m == 0 ? 1 : 2) *
                                       static_cast<double>(n == 0 ? 1 : 2) /
                                       FOUR_PI_SQUARED;
            const auto r_plus = result.spectrum.r[plus];
            const auto r_minus = result.spectrum.r[minus];
            const auto z_plus = result.spectrum.z[plus];
            const auto z_minus = result.spectrum.z[minus];
            const auto nu_plus = result.spectrum.nu[plus];
            const auto nu_minus = result.spectrum.nu[minus];
            output.rmncc[target] =
                0.5 * mode_weight * (r_plus.real() + r_minus.real());
            output.rmnss[target] =
                0.5 * mode_weight * (r_minus.real() - r_plus.real());
            output.zmnsc[target] =
                -0.5 * mode_weight * (z_plus.imag() + z_minus.imag());
            output.zmncs[target] =
                0.5 * mode_weight * (z_minus.imag() - z_plus.imag());
            output.numnsc[target] =
                -0.5 * mode_weight * (nu_plus.imag() + nu_minus.imag());
            output.numncs[target] =
                0.5 * mode_weight * (nu_minus.imag() - nu_plus.imag());
            if (output.lasym) {
                output.rmnsc[target] =
                    -0.5 * mode_weight * (r_plus.imag() + r_minus.imag());
                output.rmncs[target] =
                    0.5 * mode_weight * (r_minus.imag() - r_plus.imag());
                output.zmncc[target] =
                    0.5 * mode_weight * (z_plus.real() + z_minus.real());
                output.zmnss[target] =
                    0.5 * mode_weight * (z_minus.real() - z_plus.real());
                output.numncc[target] =
                    0.5 * mode_weight * (nu_plus.real() + nu_minus.real());
                output.numnss[target] =
                    0.5 * mode_weight * (nu_minus.real() - nu_plus.real());
            }
        }
    }

    const auto restored = restore_complex_spectrum(output, result.source_ns,
                                                   result.grid.first_surface);
    const auto parity_matches = [](const auto& actual, const auto& expected) {
        for (std::size_t index = 0; index < actual.size(); ++index) {
            const double scale = std::max(1.0, std::abs(actual[index]));
            if (std::abs(actual[index] - expected[index]) > 1.0e-10 * scale) {
                return false;
            }
        }
        return true;
    };
    if (!parity_matches(result.spectrum.r, restored.r) ||
        !parity_matches(result.spectrum.z, restored.z) ||
        !parity_matches(result.spectrum.nu, restored.nu)) {
        throw std::invalid_argument(
            result.lasym
                ? "Boozer spectra violate the real-field Hermitian contract"
                : "Boozer spectra violate the real stellarator-symmetric "
                  "parity contract");
    }
    return output;
}

MixedGridSpectrum restore_complex_spectrum(const RealBoozerSpectrum& spectrum,
                                           int source_ns,
                                           int first_surface,
                                           bool geometry,
                                           bool shift) {
    const int surfaces = source_ns - first_surface;
    validate_real_spectrum(spectrum, surfaces, geometry, shift);
    if (spectrum.mmax < 0 || spectrum.nmax < 0) {
        throw std::invalid_argument("invalid real Boozer mode bounds");
    }
    const std::size_t real_mode_count =
        checked_product(static_cast<std::size_t>(spectrum.mmax) + 1,
                        static_cast<std::size_t>(spectrum.nmax) + 1);
    const std::size_t signed_mode_count_size =
        checked_product(2 * static_cast<std::size_t>(spectrum.mmax) + 1,
                        2 * static_cast<std::size_t>(spectrum.nmax) + 1);
    if (real_mode_count > (1U << 20) || signed_mode_count_size > (1U << 22)) {
        throw std::invalid_argument("real Boozer spectrum has too many modes");
    }
    const int real_modes = static_cast<int>(real_mode_count);
    const int signed_m_count = 2 * spectrum.mmax + 1;
    const int signed_mode_count = static_cast<int>(signed_mode_count_size);
    const std::size_t expected = static_cast<std::size_t>(surfaces) *
                                 static_cast<std::size_t>(real_modes);
    if (surfaces < 1 ||
        spectrum.m.size() != static_cast<std::size_t>(real_modes) ||
        spectrum.n.size() != static_cast<std::size_t>(real_modes) ||
        (geometry && (spectrum.rmncc.size() != expected ||
                      spectrum.rmnss.size() != expected ||
                      spectrum.zmnsc.size() != expected ||
                      spectrum.zmncs.size() != expected)) ||
        (shift && (spectrum.numnsc.size() != expected ||
                   spectrum.numncs.size() != expected))) {
        throw std::invalid_argument("inconsistent real Boozer spectrum");
    }
    const auto complementary_size = spectrum.lasym ? expected : 0;
    const std::array complementary{&spectrum.rmnsc,  &spectrum.rmncs,
                                   &spectrum.zmncc,  &spectrum.zmnss,
                                   &spectrum.numncc, &spectrum.numnss};
    for (std::size_t i = 0; i < complementary.size(); ++i) {
        if (i < 4 ? !geometry : !shift) continue;
        const auto* family = complementary[i];
        if (family->size() != complementary_size) {
            throw std::invalid_argument(
                "inconsistent asymmetric Boozer spectrum");
        }
    }
    MixedGridSpectrum output;
    output.source_ns = source_ns;
    output.first_surface = first_surface;
    output.mmax = spectrum.mmax;
    output.nmax = spectrum.nmax;
    output.m.reserve(static_cast<std::size_t>(signed_mode_count));
    output.n.reserve(static_cast<std::size_t>(signed_mode_count));
    for (int n = -output.nmax; n <= output.nmax; ++n) {
        for (int m = -output.mmax; m <= output.mmax; ++m) {
            output.m.push_back(m);
            output.n.push_back(n);
        }
    }
    const std::size_t coefficient_count =
        static_cast<std::size_t>(surfaces) *
        static_cast<std::size_t>(signed_mode_count);
    if (geometry) {
        output.r.assign(coefficient_count, {});
        output.z.assign(coefficient_count, {});
    }
    if (shift) output.nu.assign(coefficient_count, {});
    const auto signed_mode = [&](int m, int n) {
        return (n + output.nmax) * signed_m_count + (m + output.mmax);
    };
    constexpr double FOUR_PI_SQUARED =
        4.0 * std::numbers::pi * std::numbers::pi;
    for (int surface = 0; surface < surfaces; ++surface) {
        for (int mode = 0; mode < real_modes; ++mode) {
            const int m = spectrum.m[static_cast<std::size_t>(mode)];
            const int n = spectrum.n[static_cast<std::size_t>(mode)];
            const int expected_m = mode % (output.mmax + 1);
            const int expected_n = mode / (output.mmax + 1);
            if (m != expected_m || n != expected_n) {
                throw std::invalid_argument("noncanonical real Boozer modes");
            }
            const std::size_t source =
                static_cast<std::size_t>(surface * real_modes + mode);
            const double integral_scale =
                FOUR_PI_SQUARED / (static_cast<double>(m == 0 ? 1 : 2) *
                                   static_cast<double>(n == 0 ? 1 : 2));
            // sin(0) contributes nothing. Keep the stored real values, but
            // do not let coincident signed modes pick up inactive terms.
            const auto coefficient = [&](const auto& family, bool sine_m,
                                         bool sine_n) {
                return family.empty() || (sine_m && m == 0) ||
                               (sine_n && n == 0)
                           ? 0.0
                           : family[source];
            };
            std::complex<double> r_plus{
                integral_scale * (coefficient(spectrum.rmncc, false, false) -
                                  coefficient(spectrum.rmnss, true, true)),
                0.0};
            std::complex<double> r_minus{
                integral_scale * (coefficient(spectrum.rmncc, false, false) +
                                  coefficient(spectrum.rmnss, true, true)),
                0.0};
            std::complex<double> z_plus{
                0.0,
                -integral_scale * (coefficient(spectrum.zmnsc, true, false) +
                                   coefficient(spectrum.zmncs, false, true))};
            std::complex<double> z_minus{
                0.0,
                integral_scale * (-coefficient(spectrum.zmnsc, true, false) +
                                  coefficient(spectrum.zmncs, false, true))};
            std::complex<double> nu_plus{
                0.0,
                -integral_scale * (coefficient(spectrum.numnsc, true, false) +
                                   coefficient(spectrum.numncs, false, true))};
            std::complex<double> nu_minus{
                0.0,
                integral_scale * (-coefficient(spectrum.numnsc, true, false) +
                                  coefficient(spectrum.numncs, false, true))};
            if (spectrum.lasym) {
                r_plus.imag(-integral_scale *
                            (coefficient(spectrum.rmnsc, true, false) +
                             coefficient(spectrum.rmncs, false, true)));
                r_minus.imag(integral_scale *
                             (-coefficient(spectrum.rmnsc, true, false) +
                              coefficient(spectrum.rmncs, false, true)));
                z_plus.real(integral_scale *
                            (coefficient(spectrum.zmncc, false, false) -
                             coefficient(spectrum.zmnss, true, true)));
                z_minus.real(integral_scale *
                             (coefficient(spectrum.zmncc, false, false) +
                              coefficient(spectrum.zmnss, true, true)));
                nu_plus.real(integral_scale *
                             (coefficient(spectrum.numncc, false, false) -
                              coefficient(spectrum.numnss, true, true)));
                nu_minus.real(integral_scale *
                              (coefficient(spectrum.numncc, false, false) +
                               coefficient(spectrum.numnss, true, true)));
            }
            const auto target = [&](int signed_m, int signed_n) {
                return static_cast<std::size_t>(
                    surface * signed_mode_count +
                    signed_mode(signed_m, signed_n));
            };
            if (geometry) {
                output.r[target(m, n)] = r_plus;
                output.r[target(m, -n)] = r_minus;
                output.r[target(-m, -n)] = std::conj(r_plus);
                output.r[target(-m, n)] = std::conj(r_minus);
                output.z[target(m, n)] = z_plus;
                output.z[target(m, -n)] = z_minus;
                output.z[target(-m, -n)] = std::conj(z_plus);
                output.z[target(-m, n)] = std::conj(z_minus);
            }
            if (shift) {
                output.nu[target(m, n)] = nu_plus;
                output.nu[target(m, -n)] = nu_minus;
                output.nu[target(-m, -n)] = std::conj(nu_plus);
                output.nu[target(-m, n)] = std::conj(nu_minus);
            }
        }
    }
    return output;
}

}  // namespace boozer_detail
}  // namespace equilibrium

namespace equilibrium {
void validate_boozer_result(const BoozerResult& result) {
    boozer_detail::validate_result(result);
}
}  // namespace equilibrium

namespace equilibrium {
Status validate_boozer_coordinates(const BoozerView& v) {
    if (v.source_ns < 2 || v.grid.first_surface < 1 ||
        v.grid.first_surface >= v.source_ns)
        return Status("invalid Boozer source radial dimensions");
    const int count = v.source_ns - v.grid.first_surface;
    if (v.nfp < 1 || v.grid.ntheta < 1 || v.grid.nzeta < 1 ||
        v.s.size() != static_cast<std::size_t>(count))
        return Status("invalid Boozer coordinate grid");
    for (std::size_t i = 0; i < v.s.size(); ++i)
        if (!std::isfinite(v.s[i]) || v.s[i] <= 0 || v.s[i] > 1 ||
            (i && v.s[i] <= v.s[i - 1]))
            return Status("Boozer radial coordinates must increase in (0,1]");
    return {};
}

Status validate_quantity(const BoozerView& v, Quantity q) {
    auto coordinates = validate_boozer_coordinates(v);
    if (!coordinates) return coordinates;
    const auto surfaces = v.s.size();
    const auto points =
        FieldGrid{static_cast<int>(surfaces), v.grid.ntheta, v.grid.nzeta}
            .size();
    if (!points) return Status("Boozer dimensions overflow");
    const auto valid = [](std::span<const double> a, std::size_t size) {
        return a.size() == size && std::ranges::all_of(a, [](double x) {
                   return std::isfinite(x);
               });
    };
    if (!v.has_quantity(q))
        return Status(std::string("missing ") + quantity_name(q));
    if (q == Quantity::GRID) return {};
    if (q == Quantity::GEOMETRY || q == Quantity::SPECTRAL_SHIFT) {
        const bool shape = q == Quantity::GEOMETRY;
        if (v.real_spectrum) {
            const auto& r = *v.real_spectrum;
            if (r.lasym != v.lasym || r.mmax < 0 || r.nmax < 0)
                return Status("invalid real spectral metadata");
            auto modes = checked_mul(static_cast<std::size_t>(r.mmax) + 1,
                                     static_cast<std::size_t>(r.nmax) + 1);
            auto count = modes ? checked_mul(*modes, surfaces) : std::nullopt;
            if (!count || *modes > (1U << 20))
                return Status("real spectrum dimensions overflow");
            const std::array families{r.rmncc,  r.rmnss,  r.zmnsc,  r.zmncs,
                                      r.numnsc, r.numncs, r.rmnsc,  r.rmncs,
                                      r.zmncc,  r.zmnss,  r.numncc, r.numnss};
            for (std::size_t i = 0; i < families.size(); ++i) {
                if ((i % 6 < 4) != shape) continue;
                if (i >= 6 && !v.lasym) {
                    if (!families[i].empty() &&
                        (!valid(families[i], *count) ||
                         !std::ranges::all_of(families[i],
                                              [](double x) { return x == 0; })))
                        return Status("asymmetric coefficients require lasym");
                } else if (!valid(families[i], *count))
                    return Status("invalid real coefficient array");
            }
            // An in-memory borrowed view may specify the canonical layout
            // implicitly. Serialized mode arrays are always supplied.
            if (!r.m.empty() || !r.n.empty()) {
                if (r.m.size() != *modes || r.n.size() != *modes)
                    return Status("invalid real mode arrays");
                for (std::size_t k = 0; k < *modes; ++k)
                    if (r.m[k] != static_cast<int>(k % (r.mmax + 1)) ||
                        r.n[k] != static_cast<int>(k / (r.mmax + 1)))
                        return Status("noncanonical real Boozer modes");
            }
        }
        if (v.has_complex_spectrum()) {
            const auto& c = v.spectrum;
            if (c.mmax < 0 || c.nmax < 0)
                return Status("invalid complex mode bounds");
            auto modes = checked_mul(2 * static_cast<std::size_t>(c.mmax) + 1,
                                     2 * static_cast<std::size_t>(c.nmax) + 1);
            auto count = modes ? checked_mul(*modes, surfaces) : std::nullopt;
            if (!count || *modes > (1U << 22) || c.m.size() != *modes ||
                c.n.size() != *modes)
                return Status("invalid complex spectrum dimensions");
            const std::array families{c.r, c.z, c.nu};
            for (std::size_t i = 0; i < families.size(); ++i) {
                if ((i < 2) != shape) continue;
                const auto values = families[i];
                if (values.size() != *count ||
                    !std::ranges::all_of(values, [](auto x) {
                        return std::isfinite(x.real()) &&
                               std::isfinite(x.imag());
                    }))
                    return Status("invalid complex coefficient array");
            }
            for (std::size_t k = 0; k < *modes; ++k)
                if (c.m[k] != static_cast<int>(k % (2 * c.mmax + 1)) - c.mmax ||
                    c.n[k] != static_cast<int>(k / (2 * c.mmax + 1)) - c.nmax)
                    return Status("noncanonical complex Boozer modes");
        }
        return {};
    }
    const auto array = [&](std::span<const double> a,
                           std::size_t count) -> Status {
        return valid(a, count) ? Status{}
                               : Status(std::string("invalid ") +
                                        quantity_name(q) + " array");
    };
    const auto orientation = [&](std::span<const double> a,
                                 std::size_t count) -> Status {
        auto status = array(a, count);
        if (!status) return status;
        for (double x : a)
            if (x == 0 || std::signbit(x) != std::signbit(a.front()))
                return Status(
                    "Jacobian orientation must be nonzero and consistent");
        return {};
    };
    if (q == Quantity::ROTATIONAL_TRANSFORM) return array(v.iota, surfaces);
    if (q == Quantity::FIELD_STRENGTH) return array(v.grid.b, *points);
    if (q == Quantity::STRENGTH_JACOBIAN) return array(v.grid.sqrtg_b, *points);
    if (q == Quantity::B2_JACOBIAN) return orientation(v.grid.b2j00, surfaces);
    if (q == Quantity::TOROIDAL_SHIFT) return array(v.grid.nu, *points);
    if (q == Quantity::SOURCE_THETA_SHIFT)
        return array(v.grid.source_theta_shift, *points);
    if (q == Quantity::PRESSURE) return array(v.profiles->pressure, surfaces);
    if (q == Quantity::SURFACE_CURRENTS) {
        if (!valid(v.profiles->current_i, surfaces) ||
            !valid(v.profiles->current_g, surfaces))
            return Status("invalid surface current arrays");
        return {};
    }
    if (q == Quantity::MAGNETIC_JACOBIAN ||
        q == Quantity::CONTRAVARIANT_FIELD || q == Quantity::COVARIANT_FIELD) {
        const auto& fields = v.magnetic_data()->get().fields;
        for (std::size_t i = 0; i < fields.size(); ++i)
            if (half_field_quantity(i) == q && !valid(fields[i], *points))
                return Status("invalid magnetic component array");
        if (q == Quantity::MAGNETIC_JACOBIAN)
            return orientation(fields[0], *points);
        return {};
    }
    if (q == Quantity::FLUX_INPUT) {
        const auto& t = *v.tracing;
        if (t.aphi.size() > (1U << 20) || !std::isfinite(t.phiedge) ||
            t.phiedge == 0 || !valid(t.aphi, t.aphi.size()))
            return Status("invalid flux input");
        try {
            profile_detail::validate_flux_type(t.pflux_type, t.aphi);
            double total = 0;
            if (t.pflux_type == "list")
                total = profile_detail::evaluate_list(t.aphi, 1., true);
            else
                for (double a : t.aphi) total += a;
            if (!std::isfinite(total) || total == 0)
                return Status("invalid toroidal flux profile");
        } catch (const std::exception& e) { return Status(e.what()); }
        return {};
    }
    if (q == Quantity::FLUX_DERIVATIVES) {
        const auto& t = *v.tracing;
        if (!valid(t.toroidal_flux_derivative, v.source_ns) ||
            !valid(t.poloidal_flux_derivative, v.source_ns))
            return Status("invalid flux derivative arrays");
        for (int i = 0; i < v.source_ns; ++i)
            if (t.toroidal_flux_derivative[i] == 0 ||
                !std::isfinite(t.poloidal_flux_derivative[i] /
                               t.toroidal_flux_derivative[i]))
                return Status("invalid signed flux derivative");
        for (std::size_t r = 0; r < surfaces; ++r)
            if (std::abs(v.s[r] -
                         full_radial_s(v.source_ns, r + v.grid.first_surface)) >
                64 * std::numeric_limits<double>::epsilon())
                return Status(
                    "flux derivatives and radial coordinates disagree");
        return {};
    }
    if (q == Quantity::AXIS_GEOMETRY) {
        if (v.source_ntor < 0) return Status("invalid axis mode count");
        const auto& t = *v.tracing;
        for (auto a : {t.axis_r_cos, t.axis_r_sin, t.axis_z_cos, t.axis_z_sin})
            if (!valid(a, static_cast<std::size_t>(v.source_ntor) + 1))
                return Status("invalid axis geometry array");
        return {};
    }
    if (q == Quantity::AXIS_FIELDS) {
        const int width = static_cast<int>(v.radial_order);
        if ((width != 2 && width != 4) || v.source_ns <= width)
            return Status("invalid axis interpolation width");
        const auto& t = *v.tracing;
        if (!valid(t.axis_half_s, width))
            return Status("invalid axis radial grid");
        for (int i = 0; i < width; ++i)
            if (t.axis_half_s[i] != half_radial_s(v.source_ns, i))
                return Status("noncanonical axis half grid");
        const auto count =
            FieldGrid{width, v.source_ntheta, v.source_nzeta}.size();
        if (!count) return Status("invalid axis field dimensions");
        for (auto a : t.axis_half_fields)
            if (!valid(a, *count)) return Status("invalid axis field array");
        auto status = orientation(t.axis_half_fields[0], *count);
        if (!status) return status;
        const auto& f = t.axis_half_fields;
        for (std::size_t i = 0; i < *count; ++i) {
            const double b2 =
                f[1][i] * f[4][i] + f[2][i] * f[5][i] + f[3][i] * f[6][i];
            if (!(b2 > 0) || !std::isfinite(b2))
                return Status("invalid axis magnetic components");
        }
        return {};
    }
    return Status("unsupported Boozer quantity");
}

std::span<const Quantity> consistency_quantities(BoozerConsistency check) {
    using Q = Quantity;
    static constexpr std::array flux{Q::FLUX_DERIVATIVES,
                                     Q::ROTATIONAL_TRANSFORM};
    static constexpr std::array currents{Q::SURFACE_CURRENTS,
                                         Q::ROTATIONAL_TRANSFORM};
    static constexpr std::array magnetic{Q::CONTRAVARIANT_FIELD,
                                         Q::COVARIANT_FIELD};
    static constexpr std::array strength{Q::FIELD_STRENGTH,
                                         Q::STRENGTH_JACOBIAN, Q::B2_JACOBIAN};
    switch (check) {
        case BoozerConsistency::FLUX_IOTA:
            return flux;
        case BoozerConsistency::SURFACE_CURRENTS:
            return currents;
        case BoozerConsistency::MAGNETIC:
            return magnetic;
        case BoozerConsistency::STRENGTH:
            return strength;
    }
    return {};
}
Status validate_consistency(const BoozerView& v, BoozerConsistency check) {
    for (auto q : consistency_quantities(check)) {
        if (!v.has_quantity(q)) return {};
        auto status = validate_quantity(v, q);
        if (!status) return status;
    }
    if (check == BoozerConsistency::FLUX_IOTA)
        for (std::size_t i = 0; i < v.iota.size(); ++i) {
            const auto r = i + v.grid.first_surface;
            const double ratio = v.tracing->poloidal_flux_derivative[r] /
                                 v.tracing->toroidal_flux_derivative[r];
            if (std::abs(v.iota[i] - ratio) >
                64 * std::numeric_limits<double>::epsilon() *
                    std::max(1., std::abs(ratio)))
                return Status("flux derivatives and iota disagree");
        }
    if (check == BoozerConsistency::SURFACE_CURRENTS)
        for (std::size_t i = 0; i < v.iota.size(); ++i) {
            const double term = v.iota[i] * v.profiles->current_i[i],
                         g = v.profiles->current_g[i];
            if (!std::isfinite(term + g) ||
                std::abs(term + g) <=
                    64 * std::numeric_limits<double>::epsilon() *
                        (std::abs(term) + std::abs(g)))
                return Status(
                    "Boozer profile G+iota*I is zero or ill-conditioned");
        }
    if (check == BoozerConsistency::MAGNETIC) {
        const auto& f = v.magnetic_data()->get().fields;
        for (std::size_t i = 0; i < f[1].size(); ++i) {
            const double b2 =
                f[1][i] * f[4][i] + f[2][i] * f[5][i] + f[3][i] * f[6][i];
            if (!(b2 > 0) || !std::isfinite(b2))
                return Status("invalid magnetic components");
        }
    }
    if (check == BoozerConsistency::STRENGTH) {
        const auto normalized = v.b2_jacobian().value();
        const auto plane = v.grid.b.size() / v.s.size();
        for (std::size_t i = 0; i < v.grid.b.size(); ++i) {
            const double expected = normalized[i / plane];
            const double actual = v.grid.b[i] * v.grid.b[i] * v.grid.sqrtg_b[i];
            if (!(v.grid.b[i] > 0) || !std::isfinite(actual) || expected == 0 ||
                std::abs(actual - expected) > 1e-9 * std::abs(expected))
                return Status(
                    "Boozer strength and Jacobian normalization disagree");
        }
    }
    return {};
}
}  // namespace equilibrium
