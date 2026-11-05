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
std::size_t checked_product(std::size_t left, std::size_t right) {
    if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left) {
        throw std::invalid_argument("Boozer result dimension overflow");
    }
    return left * right;
}

}  // namespace

namespace boozer_detail {

void validate_result(const BoozerResult& result) {
    const int surfaces = result.source_ns - result.grid.first_surface;
    if (result.spectrum.mmax < 0 || result.spectrum.nmax < 0 || surfaces < 1) {
        throw std::invalid_argument("inconsistent Boozer result dimensions");
    }
    const std::size_t mode_count =
        checked_product(2 * static_cast<std::size_t>(result.spectrum.mmax) + 1,
                        2 * static_cast<std::size_t>(result.spectrum.nmax) + 1);
    if (mode_count > (1U << 22)) {
        throw std::invalid_argument("Boozer result has too many modes");
    }
    const std::size_t points =
        checked_product(static_cast<std::size_t>(result.grid.ntheta),
                        static_cast<std::size_t>(result.grid.nzeta));
    const std::size_t real_count =
        checked_product(static_cast<std::size_t>(surfaces), points);
    const std::size_t spectral_count =
        checked_product(static_cast<std::size_t>(surfaces), mode_count);
    if (result.source_format_version < 1 || result.source_ns < 2 ||
        result.nfp < 1 || result.grid.first_surface < 1 ||
        result.grid.ntheta < 1 || result.grid.nzeta < 1 ||
        result.s.size() != static_cast<std::size_t>(surfaces) ||
        result.iota.size() != static_cast<std::size_t>(surfaces) ||
        result.grid.b2j00.size() != static_cast<std::size_t>(surfaces) ||
        result.grid.b.size() != real_count ||
        result.grid.sqrtg_b.size() != real_count ||
        result.spectrum.m.size() != mode_count ||
        result.spectrum.n.size() != mode_count ||
        result.spectrum.r.size() != spectral_count ||
        result.spectrum.z.size() != spectral_count ||
        result.spectrum.nu.size() != spectral_count) {
        throw std::invalid_argument("inconsistent Boozer result dimensions");
    }
    const auto finite_reals = [](std::span<const double> values) {
        return std::ranges::all_of(
            values, [](double value) { return std::isfinite(value); });
    };
    const auto finite_complex =
        [](std::span<const std::complex<double>> values) {
            return std::ranges::all_of(values, [](const auto value) {
                return std::isfinite(value.real()) &&
                       std::isfinite(value.imag());
            });
        };
    if (!finite_reals(result.s) || !finite_reals(result.iota) ||
        !finite_reals(result.grid.b) || !finite_reals(result.grid.sqrtg_b) ||
        !finite_reals(result.grid.b2j00) ||
        !finite_complex(result.spectrum.r) ||
        !finite_complex(result.spectrum.z) ||
        !finite_complex(result.spectrum.nu) ||
        !(result.resonance_tolerance > 0.0) ||
        !std::isfinite(result.resonance_tolerance)) {
        throw std::invalid_argument("Boozer result contains non-finite data");
    }
    if (result.profiles) {
        if (!result.tracing || !result.tracing->magnetic)
            throw std::invalid_argument(
                "Boozer surface profiles require v6 tracing and fields");
        const auto& profiles = *result.profiles;
        for (const auto* values :
             {&profiles.pressure, &profiles.current_i, &profiles.current_g}) {
            if (values->size() != result.s.size() || !finite_reals(*values))
                throw std::invalid_argument(
                    "invalid Boozer surface profile array");
        }
        for (std::size_t j = 0; j < result.s.size(); ++j) {
            const double iota_i = result.iota[j] * profiles.current_i[j];
            const double g = profiles.current_g[j];
            const double normalization = g + iota_i;
            const double scale = std::abs(g) + std::abs(iota_i);
            if (!std::isfinite(normalization) ||
                std::abs(normalization) <=
                    64.0 * std::numeric_limits<double>::epsilon() * scale)
                throw std::invalid_argument(
                    "Boozer profile G+iota*I is zero or ill-conditioned");
        }
    }
    if (result.tracing) {
        const auto& tracing = *result.tracing;
        const int width = static_cast<int>(result.radial_order);
        if ((width != 2 && width != 4) || result.source_ns <= width ||
            result.source_ntheta < 1 || result.source_nzeta < 1 ||
            result.source_ntor < 0 || tracing.aphi.empty() ||
            tracing.aphi.size() > (1U << 20) ||
            !std::isfinite(tracing.phiedge) || tracing.phiedge == 0.0) {
            throw std::invalid_argument("invalid Boozer tracing metadata");
        }
        const auto require = [&](std::span<const double> values,
                                 std::size_t count) {
            if (values.size() != count || !finite_reals(values))
                throw std::invalid_argument("invalid Boozer tracing array");
        };
        const auto ns = static_cast<std::size_t>(result.source_ns);
        const auto axis_modes =
            static_cast<std::size_t>(result.source_ntor) + 1;
        require(tracing.aphi, tracing.aphi.size());
        double profile_sum = 0.0;
        profile_detail::validate_flux_type(tracing.pflux_type, tracing.aphi);
        if (tracing.pflux_type == "list" && !result.profiles)
            throw std::invalid_argument(
                "list flux export requires surface profiles");
        if (tracing.pflux_type == "list")
            profile_sum =
                profile_detail::evaluate_list(tracing.aphi, 1.0, true);
        else
            for (double value : tracing.aphi) profile_sum += value;
        if (!std::isfinite(profile_sum) || profile_sum == 0.0)
            throw std::invalid_argument("invalid Boozer toroidal flux profile");
        require(tracing.toroidal_flux_derivative, ns);
        require(tracing.poloidal_flux_derivative, ns);
        require(tracing.axis_r_cos, axis_modes);
        require(tracing.axis_r_sin, axis_modes);
        require(tracing.axis_z_cos, axis_modes);
        require(tracing.axis_z_sin, axis_modes);
        require(tracing.axis_half_s, static_cast<std::size_t>(width));
        for (std::size_t j = 0; j < ns; ++j) {
            const double phip = tracing.toroidal_flux_derivative[j];
            const double iota = tracing.poloidal_flux_derivative[j] / phip;
            if (phip == 0.0 || !std::isfinite(iota))
                throw std::invalid_argument(
                    "invalid signed Boozer flux derivative");
            if (j >= static_cast<std::size_t>(result.grid.first_surface)) {
                const auto index =
                    j - static_cast<std::size_t>(result.grid.first_surface);
                const double tolerance =
                    64.0 * std::numeric_limits<double>::epsilon();
                if (std::abs(result.iota[index] - iota) >
                        tolerance * std::max(1.0, std::abs(iota)) ||
                    std::abs(result.s[index] -
                             static_cast<double>(j) /
                                 static_cast<double>(ns - 1)) > tolerance)
                    throw std::invalid_argument(
                        "inconsistent Boozer tracing profiles");
            }
        }
        for (std::size_t j = 0; j < tracing.axis_half_s.size(); ++j) {
            const double expected =
                (static_cast<double>(j) + 0.5) / static_cast<double>(ns - 1);
            if (tracing.axis_half_s[j] != expected)
                throw std::invalid_argument(
                    "noncanonical Boozer axis half grid");
        }
        const auto half_count = checked_product(
            tracing.axis_half_s.size(),
            checked_product(static_cast<std::size_t>(result.source_ntheta),
                            static_cast<std::size_t>(result.source_nzeta)));
        for (const auto& field : tracing.axis_half_fields)
            require(field, half_count);
        if (tracing.magnetic) {
            const auto& magnetic = *tracing.magnetic;
            for (const auto& field : magnetic.fields)
                require(field, real_count);
            for (std::size_t j = 0; j < real_count; ++j) {
                const auto& f = magnetic.fields;
                const double b2 = f[EquilibriumSnapshot::BSUPS][j] *
                                      f[EquilibriumSnapshot::BSUBS][j] +
                                  f[EquilibriumSnapshot::BSUPU][j] *
                                      f[EquilibriumSnapshot::BSUBU][j] +
                                  f[EquilibriumSnapshot::BSUPV][j] *
                                      f[EquilibriumSnapshot::BSUBV][j];
                if (!(b2 > 0.0) || !std::isfinite(b2) ||
                    f[EquilibriumSnapshot::SQRTG][j] == 0.0)
                    throw std::invalid_argument(
                        "invalid Boozer magnetic field");
            }
            require(result.grid.source_theta_shift, real_count);
            require(result.grid.nu, real_count);
        }
    }
    std::size_t mode = 0;
    for (int n = -result.spectrum.nmax; n <= result.spectrum.nmax; ++n) {
        for (int m = -result.spectrum.mmax; m <= result.spectrum.mmax; ++m) {
            if (result.spectrum.m[mode] != m || result.spectrum.n[mode] != n) {
                throw std::invalid_argument(
                    "Boozer result has noncanonical mode ordering");
            }
            ++mode;
        }
    }
}

RealBoozerSpectrum make_real_spectrum(const BoozerResult& result) {
    validate_result(result);
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
                                           int first_surface) {
    const int surfaces = source_ns - first_surface;
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
        spectrum.rmncc.size() != expected ||
        spectrum.rmnss.size() != expected ||
        spectrum.zmnsc.size() != expected ||
        spectrum.zmncs.size() != expected ||
        spectrum.numnsc.size() != expected ||
        spectrum.numncs.size() != expected) {
        throw std::invalid_argument("inconsistent real Boozer spectrum");
    }
    const auto complementary_size = spectrum.lasym ? expected : 0;
    for (const auto* family :
         {&spectrum.rmnsc, &spectrum.rmncs, &spectrum.zmncc, &spectrum.zmnss,
          &spectrum.numncc, &spectrum.numnss}) {
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
    output.r.assign(coefficient_count, {});
    output.z.assign(coefficient_count, {});
    output.nu.assign(coefficient_count, {});
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
            std::complex<double> r_plus{
                integral_scale *
                    (spectrum.rmncc[source] - spectrum.rmnss[source]),
                0.0};
            std::complex<double> r_minus{
                integral_scale *
                    (spectrum.rmncc[source] + spectrum.rmnss[source]),
                0.0};
            std::complex<double> z_plus{
                0.0, -integral_scale *
                         (spectrum.zmnsc[source] + spectrum.zmncs[source])};
            std::complex<double> z_minus{
                0.0, integral_scale *
                         (-spectrum.zmnsc[source] + spectrum.zmncs[source])};
            std::complex<double> nu_plus{
                0.0, -integral_scale *
                         (spectrum.numnsc[source] + spectrum.numncs[source])};
            std::complex<double> nu_minus{
                0.0, integral_scale *
                         (-spectrum.numnsc[source] + spectrum.numncs[source])};
            if (spectrum.lasym) {
                r_plus.imag(-integral_scale *
                            (spectrum.rmnsc[source] + spectrum.rmncs[source]));
                r_minus.imag(integral_scale * (-spectrum.rmnsc[source] +
                                               spectrum.rmncs[source]));
                z_plus.real(integral_scale *
                            (spectrum.zmncc[source] - spectrum.zmnss[source]));
                z_minus.real(integral_scale *
                             (spectrum.zmncc[source] + spectrum.zmnss[source]));
                nu_plus.real(integral_scale * (spectrum.numncc[source] -
                                               spectrum.numnss[source]));
                nu_minus.real(integral_scale * (spectrum.numncc[source] +
                                                spectrum.numnss[source]));
            }
            const auto target = [&](int signed_m, int signed_n) {
                return static_cast<std::size_t>(
                    surface * signed_mode_count +
                    signed_mode(signed_m, signed_n));
            };
            output.r[target(m, n)] = r_plus;
            output.r[target(m, -n)] = r_minus;
            output.r[target(-m, -n)] = std::conj(r_plus);
            output.r[target(-m, n)] = std::conj(r_minus);
            output.z[target(m, n)] = z_plus;
            output.z[target(m, -n)] = z_minus;
            output.z[target(-m, -n)] = std::conj(z_plus);
            output.z[target(-m, n)] = std::conj(z_minus);
            output.nu[target(m, n)] = nu_plus;
            output.nu[target(m, -n)] = nu_minus;
            output.nu[target(-m, -n)] = std::conj(nu_plus);
            output.nu[target(-m, n)] = std::conj(nu_minus);
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
