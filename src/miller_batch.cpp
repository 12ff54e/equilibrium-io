#include "equilibrium/miller_batch.hpp"

#include "gfile_boundary.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <numbers>
#include <stdexcept>

namespace equilibrium {
namespace {
constexpr double PERIOD = 2 * std::numbers::pi;

std::vector<double> samples(std::optional<MillerRange> range,
                            double measured,
                            std::size_t limit) {
    if (!range) return {measured};
    const auto [first, last, step] = *range;
    if (!std::isfinite(first) || !std::isfinite(last) || !std::isfinite(step) ||
        step == 0 || (first != last && ((last > first) != (step > 0))))
        throw std::invalid_argument(
            "Miller range needs finite endpoints and a nonzero step toward the "
            "end");
    if (first == last) return {first};
    const long double steps = (static_cast<long double>(last) - first) / step;
    const long double nearest = std::round(steps);
    // Subtracting close, rounded decimal endpoints can cost more than one
    // ulp of the interval count (e.g. 6.0:6.1:0.1).
    const long double scale =
        std::max(std::abs(static_cast<long double>(first)),
                 std::abs(static_cast<long double>(last)));
    const long double rounding = 4 * std::numeric_limits<double>::epsilon() *
                                 (steps + scale / std::abs(step));
    const bool reaches_end =
        nearest >= 1 && std::abs(steps - nearest) <= std::min(rounding, 1.e-8L);
    const long double intervals = reaches_end ? nearest : std::floor(steps);
    if (!std::isfinite(intervals) ||
        intervals >= static_cast<long double>(limit))
        throw std::invalid_argument("Miller batch exceeds max_cases");
    const auto count = static_cast<std::size_t>(intervals) + 1;
    std::vector<double> values;
    values.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        const double value = reaches_end && i == count - 1
                                 ? last
                                 : static_cast<double>(std::fma(
                                       static_cast<long double>(i),
                                       static_cast<long double>(step),
                                       static_cast<long double>(first)));
        if (i && ((step > 0 && value <= values.back()) ||
                  (step < 0 && value >= values.back())))
            throw std::invalid_argument(
                "Miller range step is below double resolution");
        values.push_back(value);
    }
    return values;
}

struct MillerBoundary {
    double r0, a, z0, ku, kl, du, dl;

    LcfsPoint operator()(double t) const {
        // A single sinusoidal Z in the Fourier angle gives both halves the
        // same tangent speed at each join. Recover the Miller angle on each
        // arc from Z. The locus is exactly the two specified Miller arcs.
        const double height = std::lerp(-kl, ku, 0.5 * (1 + std::sin(t)));
        const bool upper = height >= 0;
        const double s = std::clamp(height / (upper ? ku : kl), -1.0, 1.0);
        const double c =
            std::copysign(std::sqrt((1 - s) * (1 + s)), std::cos(t));
        const double shift = std::asin(upper ? du : dl) * s;
        return {r0 + a * (c * std::cos(shift) - s * std::sin(shift)),
                z0 + a * height};
    }
};

LcfsShape check_fit(const MillerBoundary& target,
                    const InputParams& p,
                    std::size_t count,
                    double tolerance) {
    std::vector<LcfsPoint> polygon(count);
    for (std::size_t k = 0; k < count; ++k) {
        const double t =
            PERIOD * static_cast<double>(k) / static_cast<double>(count);
        auto& point = polygon[k];
        double dr = 0, dz = 0;
        for (int m = 0; m < p.mpol; ++m) {
            const double c = std::cos(m * t), s = std::sin(m * t);
            const double rs = p.lasym ? p.rbsc[m] : 0;
            const double zc = p.lasym ? p.zbcc[m] : 0;
            point.r += p.rbcc[m] * c + rs * s;
            point.z += zc * c + p.zbsc[m] * s;
            dr += m * (rs * c - p.rbcc[m] * s);
            dz += m * (p.zbsc[m] * c - zc * s);
        }
        const auto exact = target(t);
        if (!(point.r > 0 &&
              (point.r - target.r0) * dz - (point.z - target.z0) * dr > 0 &&
              std::hypot(point.r - exact.r, point.z - exact.z) <=
                  tolerance * target.a))
            throw std::domain_error(
                "Miller Fourier fit fails the contour or orientation check");
    }
    const auto shape = calculate_lcfs_shape(polygon, target.z0);
    const double error =
        std::max({std::abs(shape.r_center - target.r0) / target.a,
                  std::abs(shape.minor_radius - target.a) / target.a,
                  std::abs(shape.elongation_upper - target.ku),
                  std::abs(shape.elongation_lower - target.kl),
                  std::abs(shape.triangularity_upper - target.du),
                  std::abs(shape.triangularity_lower - target.dl)});
    if (error > tolerance)
        throw std::domain_error(
            "Miller Fourier fit misses requested LCFS statistics; increase "
            "mpol or adjust tolerance");
    return shape;
}

std::string filename(std::size_t i) {
    return "case_" + std::to_string(i) + ".json";
}
}  // namespace

std::size_t write_miller_batch(const std::filesystem::path& gfile,
                               const std::filesystem::path& output_directory,
                               const MillerScan& scan,
                               const GfileInputSettings& settings) {
    if (settings.psi_ratio != 1.0)
        throw std::invalid_argument(
            "Miller batches use the full G-file LCFS and profiles "
            "(psi_ratio=1)");
    if (scan.max_cases == 0)
        throw std::invalid_argument("max_cases must be positive");
    const auto source = gfile_detail::read_profiles(gfile, settings);
    const auto base = calculate_lcfs_shape(source.lcfs);
    const std::array ranges{scan.r0, scan.kappa_upper, scan.kappa_lower,
                            scan.delta_upper, scan.delta_lower};
    const std::array defaults{base.r_center, base.elongation_upper,
                              base.elongation_lower, base.triangularity_upper,
                              base.triangularity_lower};
    std::array<std::vector<double>, 5> grid;
    std::size_t total = 1;
    for (std::size_t i = 0; i < grid.size(); ++i) {
        grid[i] = samples(ranges[i], defaults[i], scan.max_cases / total);
        total *= grid[i].size();
    }
    for (double r : grid[0])
        if (!(r > base.minor_radius))
            throw std::domain_error("Miller R0 must exceed minor radius");
    for (auto i : {1, 2})
        for (double k : grid[i])
            if (!(k > 0))
                throw std::domain_error("Miller elongation must be positive");
    for (auto i : {3, 4})
        for (double d : grid[i])
            if (!(std::abs(d) < 1))
                throw std::domain_error(
                    "Miller triangularity must lie strictly between -1 and 1");

    // Round up so the DFT samples always include the top and bottom.
    const std::size_t n = ((settings.boundary_points + 3) / 4) * 4;
    std::vector<LcfsPoint> points(n);
    const auto index_path = output_directory / "index.csv";
    const auto partial = output_directory / "index.csv.pending";
    if (!std::filesystem::create_directory(output_directory))
        throw std::runtime_error("Miller output directory already exists: " +
                                 output_directory.string());
    std::size_t attempted = 0;
    try {
        std::ofstream index;
        index.exceptions(std::ios::badbit | std::ios::failbit);
        index.open(partial);
        index.imbue(std::locale::classic());
        index << std::setprecision(17)
              << "file,r0,kappa_upper,kappa_lower,delta_upper,delta_lower,"
                 "minor_radius,z_reference,"
                 "fitted_r0,fitted_minor_radius,fitted_kappa_upper,fitted_"
                 "kappa_lower,fitted_delta_upper,fitted_delta_lower\n";
        for (std::size_t i = 0; i < total; ++i) {
            std::array<double, 5> values;
            auto pos = i;
            for (std::size_t d = grid.size(); d-- > 0;) {
                values[d] = grid[d][pos % grid[d].size()];
                pos /= grid[d].size();
            }
            const MillerBoundary target{
                values[0], base.minor_radius, base.z_reference, values[1],
                values[2], values[3],         values[4]};
            for (std::size_t k = 0; k < n; ++k)
                points[k] = target(PERIOD * static_cast<double>(k) /
                                   static_cast<double>(n));
            GfileInput result;
            result.input = source.input;
            try {
                gfile_detail::fit_boundary(result.input, points,
                                           {target.r0, target.z0}, target.a,
                                           settings.boundary_tolerance);
                const auto measured = check_fit(
                    target, result.input, std::max<std::size_t>(4096, 2 * n),
                    settings.boundary_tolerance);
                const auto name = filename(i);
                attempted = i + 1;
                const auto written = write_gfile_input(
                    (output_directory / name).string(), result);
                if (!written) throw std::runtime_error(written.error());
                index << name;
                for (double x : values) index << ',' << x;
                index << ',' << target.a << ',' << target.z0 << ','
                      << measured.r_center << ',' << measured.minor_radius
                      << ',' << measured.elongation_upper << ','
                      << measured.elongation_lower << ','
                      << measured.triangularity_upper << ','
                      << measured.triangularity_lower << '\n';
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "Miller case " + std::to_string(i) +
                    " (R0=" + std::to_string(target.r0) +
                    ", kappa_upper/lower=" + std::to_string(target.ku) + "/" +
                    std::to_string(target.kl) +
                    ", delta_upper/lower=" + std::to_string(target.du) + "/" +
                    std::to_string(target.dl) + ", tolerance=" +
                    std::to_string(settings.boundary_tolerance) +
                    "): " + e.what());
            }
        }
        index.close();
        std::filesystem::rename(partial, index_path);
    } catch (...) {
        std::error_code error;
        std::filesystem::remove(partial, error);
        for (std::size_t i = 0; i < attempted; ++i)
            std::filesystem::remove(output_directory / filename(i), error);
        std::filesystem::remove(output_directory, error);
        throw;
    }
    return total;
}
}  // namespace equilibrium
