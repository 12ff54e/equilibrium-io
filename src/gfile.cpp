#include "equilibrium/gfile.hpp"

#include "equilibrium/record_validation.hpp"
#include "equilibrium/writer_helpers.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <memory>
#include <numbers>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string_view>

#include <Contour.h>
#include <GFileRawData.h>

namespace equilibrium {
namespace {
constexpr double PERIOD = 2.0 * std::numbers::pi;
constexpr std::size_t MAX_FILE_BYTES = 128U << 20;
struct CloseFile {
    void operator()(std::FILE* file) const { std::fclose(file); }
};

void require(bool condition, const char* message) {
    if (!condition) throw std::invalid_argument(message);
}

bool same_point(LcfsPoint a, LcfsPoint b) {
    return a.r == b.r && a.z == b.z;
}

std::vector<LcfsPoint> clean_polygon(std::span<const LcfsPoint> points) {
    std::vector<LcfsPoint> clean;
    for (auto point : points) {
        require(
            std::isfinite(point.r) && point.r > 0.0 && std::isfinite(point.z),
            "G-file LCFS needs finite coordinates and positive R");
        if (clean.empty() || !same_point(point, clean.back()))
            clean.push_back(point);
    }
    if (clean.size() > 1 && same_point(clean.front(), clean.back()))
        clean.pop_back();
    require(clean.size() >= 3,
            "G-file LCFS needs at least three distinct points");
    return clean;
}

int header_integer(std::string_view text) {
    const auto first = text.find_first_not_of(' ');
    require(first != std::string_view::npos, "missing G-file grid dimension");
    text.remove_prefix(first);
    const auto last = text.find_last_not_of(' ');
    text = text.substr(0, last + 1);
    int value = 0;
    const auto [end, error] =
        std::from_chars(text.data(), text.data() + text.size(), value);
    require(error == std::errc{} && end == text.data() + text.size(),
            "invalid G-file grid dimension");
    return value;
}

// libmeq owns format interpretation and COCOS inference. Check the bounded
// numeric payload first: its reader assumes valid counts and complete arrays.
GFileRawData<double> load_gfile(std::istream& input, std::optional<int> cocos) {
    std::optional<Cocos> convention;
    if (cocos) convention.emplace(*cocos);
    std::string text;
    std::array<char, 8192> chunk{};
    while (
        input.read(chunk.data(), static_cast<std::streamsize>(chunk.size())) ||
        input.gcount() != 0) {
        require(text.size() + static_cast<std::size_t>(input.gcount()) <=
                    MAX_FILE_BYTES,
                "G-file exceeds the 128 MiB input limit");
        text.append(chunk.data(), static_cast<std::size_t>(input.gcount()));
    }
    require(input.eof() && !input.bad(), "cannot read G-file stream");
    const auto newline = text.find('\n');
    require(newline != std::string::npos && newline >= 60 &&
                text.find('\0') == std::string::npos,
            "G-file needs a 48-character header and three I4 fields");
    const int nw = header_integer(std::string_view(text).substr(52, 4));
    const int nh = header_integer(std::string_view(text).substr(56, 4));
    require(nw >= 2 && nh >= 2 && nw <= 4096 && nh <= 4096 &&
                static_cast<std::size_t>(nw) * static_cast<std::size_t>(nh) <=
                    (1U << 22),
            "G-file grid dimensions must be 2-4096 with at most 2^22 points");
    for (std::size_t k = newline + 1; k < text.size(); ++k)
        if (text[k] == 'D' || text[k] == 'd') text[k] = 'E';
    std::istringstream preflight(text.substr(newline + 1));
    preflight.imbue(std::locale::classic());
    const auto number = [](std::istream& stream) {
        double value = 0.0;
        require(static_cast<bool>(stream >> value) && std::isfinite(value),
                "G-file contains truncated or nonfinite numeric data");
        return value;
    };
    std::array<double, 20> scalars{};
    for (std::size_t row = 0; row < 4; ++row) {
        std::string line;
        require(static_cast<bool>(std::getline(preflight, line)),
                "G-file scalar header is truncated");
        std::istringstream values(line);
        values.imbue(std::locale::classic());
        for (std::size_t column = 0; column < 5; ++column)
            scalars[5 * row + column] = number(values);
        values >> std::ws;
        require(values.eof(),
                "G-file scalar header must have five values per row");
    }
    const LcfsPoint axis{scalars[5], scalars[6]};
    require(scalars[0] > 0.0 && scalars[1] > 0.0 && scalars[3] > 0.0 &&
                axis.r > 0.0 && std::isfinite(scalars[3] + scalars[0]) &&
                std::isfinite(scalars[4] - .5 * scalars[1]) &&
                std::isfinite(scalars[4] + .5 * scalars[1]),
            "invalid G-file grid extent or magnetic axis");
    const std::size_t array_count =
        static_cast<std::size_t>(nw) * (static_cast<std::size_t>(nh) + 5);
    for (std::size_t k = 0; k < array_count; ++k) (void)number(preflight);
    long long boundary_count = 0, limiter_count = 0;
    require(static_cast<bool>(preflight >> boundary_count >> limiter_count) &&
                boundary_count >= 3 && boundary_count <= 100000 &&
                limiter_count >= 0 && limiter_count <= 100000,
            "invalid G-file boundary or limiter count");
    std::vector<LcfsPoint> boundary;
    boundary.reserve(static_cast<std::size_t>(boundary_count));
    for (long long k = 0; k < boundary_count; ++k)
        boundary.push_back({number(preflight), number(preflight)});
    boundary = clean_polygon(boundary);
    // Avoid the reader's unchecked access after deduplicating polar angles.
    std::vector<double> angles;
    for (auto p : boundary) {
        require(p.r != axis.r || p.z != axis.z,
                "G-file boundary includes the magnetic axis");
        // Match libmeq's angle normalization exactly: adding 2*pi can merge
        // distinct tiny negative atan2 values through rounding.
        angles.push_back(util::arctan(Vec<2, double>{p.r, p.z} -
                                      Vec<2, double>{axis.r, axis.z}));
    }
    std::sort(angles.begin(), angles.end());
    require(std::unique(angles.begin(), angles.end()) - angles.begin() >= 3,
            "degenerate G-file boundary about the magnetic axis");
    for (long long k = 0; k < 2 * limiter_count; ++k) (void)number(preflight);

    std::istringstream stream(text);
    stream.imbue(std::locale::classic());
    GFileRawData<double> data;
    data.read(stream, convention);
    require(data.is_complete(), "libmeq could not read G-file data");
    // Preserve the supplied polygon order; libmeq rearranges it for Boozer.
    data.boundary.clear();
    for (auto p : boundary) data.boundary.emplace_back(p.r, p.z);
    data.boundary_num = static_cast<unsigned>(boundary.size());
    return data;
}

GfileLcfs lcfs_of(const GFileRawData<double>& data) {
    GfileLcfs result;
    result.magnetic_axis = {data.magnetic_axis.x(), data.magnetic_axis.y()};
    result.cocos = data.cocos.index();
    for (auto p : data.boundary) result.points.push_back({p.x(), p.y()});
    return result;
}

std::ifstream open_gfile(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream)
        throw std::runtime_error("cannot open G-file: " + path.string());
    return stream;
}

void validate_settings(const GfileInputSettings& settings) {
    require(settings.mpol >= 2 && settings.mpol <= 16,
            "G-file conversion mpol must be 2-16");
    require(settings.ntheta == 0 ||
                (settings.ntheta >= 2 * settings.mpol + 6 &&
                 settings.ntheta <= 256 && settings.ntheta % 2 == 0),
            "G-file ntheta must be even, at least 2*mpol+6, and at most 256");
    require(settings.boundary_points >=
                    static_cast<std::size_t>(4 * settings.mpol) &&
                settings.boundary_points <= 16384 &&
                settings.profile_points >= 2 &&
                settings.profile_points <= 65536,
            "invalid G-file boundary/profile sampling dimensions");
    require(std::isfinite(settings.psi_ratio) && settings.psi_ratio > 0.0 &&
                settings.psi_ratio <= 1.0,
            "G-file psi_ratio must be in (0,1]");
    require(std::isfinite(settings.boundary_tolerance) &&
                settings.boundary_tolerance > 0.0 &&
                std::isfinite(settings.poloidal_flux_tolerance) &&
                settings.poloidal_flux_tolerance > 0.0,
            "G-file conversion tolerances must be positive and finite");
    require(settings.stage.ns >= 3 && settings.stage.ns <= 512 &&
                settings.stage.max_iter >= 1 &&
                std::isfinite(settings.stage.ftol) &&
                settings.stage.ftol >= 1.e-16,
            "G-file stage needs ns in [3,512], max_iter>=1, and ftol>=1e-16");
}

// Piecewise linear q on normalized source poloidal flux, integrated exactly.
class FluxProfile {
   public:
    FluxProfile(std::vector<double> q, double limit) : q_(std::move(q)) {
        const double sign = std::copysign(1.0, q_.front());
        for (double value : q_)
            require(value * sign > 0.0 && std::isfinite(value),
                    "G-file q must stay nonzero with one sign");
        integral_.resize(q_.size());
        const double dx = 1.0 / static_cast<double>(q_.size() - 1);
        for (std::size_t k = 1; k < q_.size(); ++k) {
            integral_[k] = integral_[k - 1] + 0.5 * (q_[k - 1] + q_[k]) * dx;
            require(std::isfinite(integral_[k]) &&
                        sign * integral_[k] > sign * integral_[k - 1],
                    "G-file toroidal-flux map overflows or loses monotonicity");
        }
        limit_ = limit;
        edge_ = integral(limit_);
    }
    static double interpolate(const std::vector<double>& values, double x) {
        const double cell = x * static_cast<double>(values.size() - 1);
        const auto k =
            std::min(static_cast<std::size_t>(cell), values.size() - 2);
        return std::lerp(values[k], values[k + 1],
                         cell - static_cast<double>(k));
    }
    double q(double x) const { return interpolate(q_, x); }
    double integral(double x) const {
        const double cell = x * static_cast<double>(q_.size() - 1);
        const auto k = std::min(static_cast<std::size_t>(cell), q_.size() - 2);
        const double dx =
            x - static_cast<double>(k) / static_cast<double>(q_.size() - 1);
        return integral_[k] + 0.5 * (q_[k] + q(x)) * dx;
    }
    double edge() const { return edge_; }
    double invert(double toroidal_fraction) const {
        if (toroidal_fraction == 0.0) return 0.0;
        if (toroidal_fraction == 1.0) return limit_;
        const double target = toroidal_fraction * edge_;
        return util::find_root(
            [&](double x) { return integral(x) - target; }, 0.0, limit_,
            [](double a, double b) {
                return b - a < 8.0 * std::numeric_limits<double>::epsilon();
            });
    }

   private:
    std::vector<double> q_, integral_;
    double limit_ = 1.0, edge_ = 0.0;
};

std::vector<LcfsPoint> sample_boundary(GFileRawData<double>& data,
                                       const GfileInputSettings& settings) {
    const auto axis = data.magnetic_axis;
    const double left = data.r_left, right = left + data.dim.x();
    const double bottom = data.z_mid - 0.5 * data.dim.y();
    const double top = data.z_mid + 0.5 * data.dim.y();
    const auto inside = [&](auto point) {
        return point.x() >= left && point.x() <= right && point.y() >= bottom &&
               point.y() <= top;
    };
    require(inside(axis), "G-file magnetic axis is outside the flux grid");
    std::vector<Vec<2, double>> contour = data.boundary;
    double winding = 0.0;
    for (std::size_t k = 0; k < contour.size(); ++k) {
        require(inside(contour[k]), "G-file LCFS is outside the flux grid");
        const auto a = contour[k] - axis;
        const auto b = contour[(k + 1) % contour.size()] - axis;
        const double step = std::atan2(a.x() * b.y() - a.y() * b.x(),
                                       a.x() * b.x() + a.y() * b.y());
        require(std::abs(step) > 1.e-12 && std::abs(step) < std::numbers::pi &&
                    (k == 0 || step * winding > 0.0),
                "G-file conversion requires an ordered star-shaped LCFS about "
                "the axis");
        winding += step;
    }
    require(std::abs(std::abs(winding) - PERIOD) < 1.e-8,
            "G-file LCFS must wind once around the magnetic axis");
    if (settings.psi_ratio < 1.0) {
        require(data.nw >= 6 && data.nh >= 6,
                "G-file contour cropping requires at least a 6x6 flux grid");
        const double delta = data.flux_LCFS - data.flux_magnetic_axis;
        auto normalized_flux = data.flux;
        for (unsigned r = 0; r < data.nw; ++r)
            for (unsigned z = 0; z < data.nh; ++z)
                normalized_flux(r, z) =
                    (normalized_flux(r, z) - data.flux_magnetic_axis) / delta;
        const intp::InterpolationFunction<double, 2, 5> flux(
            normalized_flux, std::pair{left, right}, std::pair{bottom, top});
        require(
            std::abs(flux(axis)) <= 1.e-4 && flux(axis) < settings.psi_ratio,
            "interpolated G-file axis flux disagrees with the header or crop");
        for (auto point : contour) {
            require(flux(point) > settings.psi_ratio,
                    "requested G-file crop is not bracketed by every LCFS ray");
        }
        const Contour<double> cropped(settings.psi_ratio, flux, data);
        for (std::size_t k = 0; k < contour.size(); ++k) {
            contour[k] = cropped[k];
            require(std::abs(flux(contour[k]) - settings.psi_ratio) <= 1.e-9,
                    "G-file cropped contour failed its flux residual check");
        }
    }
    std::vector<std::pair<double, double>> polar;
    for (auto point : contour) {
        const auto relative = point - axis;
        double angle = std::atan2(relative.y(), relative.x());
        if (angle < 0.0) angle += PERIOD;
        polar.emplace_back(angle, std::hypot(relative.x(), relative.y()));
    }
    std::sort(polar.begin(), polar.end());
    require(polar.size() >= 4,
            "G-file boundary interpolation needs four points");
    std::vector<double> angle, radius;
    for (auto [t, r] : polar) {
        angle.push_back(t);
        radius.push_back(r);
    }
    angle.push_back(angle.front() + PERIOD);
    const intp::InterpolationFunction<double, 1, 3> radial(
        true, std::pair{radius.begin(), radius.end()},
        std::pair{angle.begin(), angle.end()});
    std::vector<LcfsPoint> sampled(settings.boundary_points);
    for (std::size_t k = 0; k < sampled.size(); ++k) {
        const double t = PERIOD * static_cast<double>(k) /
                         static_cast<double>(sampled.size());
        const double r = radial(t < angle.front() ? t + PERIOD : t);
        require(std::isfinite(r) && r > 0.0,
                "G-file boundary interpolation has a nonpositive radius");
        sampled[k] = {axis.x() + r * std::cos(t), axis.y() + r * std::sin(t)};
        require(sampled[k].r > 0.0, "interpolated G-file LCFS crosses R=0");
    }
    return sampled;
}

void fit_boundary(InputParams& input,
                  std::span<const LcfsPoint> sampled,
                  GfileInput& result,
                  const GfileInputSettings& settings) {
    std::vector<double> rc(input.mpol), rs(input.mpol), zc(input.mpol),
        zs(input.mpol);
    for (int m = 0; m < input.mpol; ++m) {
        for (std::size_t k = 0; k < sampled.size(); ++k) {
            const double t = PERIOD * static_cast<double>(k) /
                             static_cast<double>(sampled.size());
            const double c = std::cos(m * t), s = std::sin(m * t);
            rc[m] += sampled[k].r * c;
            rs[m] += sampled[k].r * s;
            zc[m] += sampled[k].z * c;
            zs[m] += sampled[k].z * s;
        }
        const double factor =
            (m == 0 ? 1.0 : 2.0) / static_cast<double>(sampled.size());
        rc[m] *= factor;
        rs[m] *= factor;
        zc[m] *= factor;
        zs[m] *= factor;
    }
    const double a = calculate_lcfs_shape(result.lcfs).minor_radius;
    const double symmetry_tolerance =
        1.e-12 * std::max(a, std::abs(result.lcfs.magnetic_axis.r));
    input.lasym = std::abs(result.lcfs.magnetic_axis.z) > symmetry_tolerance;
    for (int m = 0; m < input.mpol; ++m)
        input.lasym = input.lasym || std::abs(rs[m]) > symmetry_tolerance ||
                      std::abs(zc[m]) > symmetry_tolerance;
    if (!input.lasym) {
        std::fill(rs.begin(), rs.end(), 0.0);
        std::fill(zc.begin(), zc.end(), 0.0);
    }
    for (std::size_t k = 0; k < sampled.size(); ++k) {
        const double t = PERIOD * static_cast<double>(k) /
                         static_cast<double>(sampled.size());
        double r = 0.0, z = 0.0, dr = 0.0, dz = 0.0;
        for (int m = 0; m < input.mpol; ++m) {
            const double c = std::cos(m * t), s = std::sin(m * t);
            r += rc[m] * c + rs[m] * s;
            z += zc[m] * c + zs[m] * s;
            dr += m * (-rc[m] * s + rs[m] * c);
            dz += m * (-zc[m] * s + zs[m] * c);
        }
        require(r > 0.0 && (r - result.lcfs.magnetic_axis.r) * dz -
                                   (z - result.lcfs.magnetic_axis.z) * dr >
                               0.0,
                "G-file Fourier boundary folds or crosses R=0");
        result.boundary_max_error =
            std::max(result.boundary_max_error,
                     std::hypot(r - sampled[k].r, z - sampled[k].z));
    }
    require(result.boundary_max_error <= settings.boundary_tolerance * a,
            "G-file Fourier boundary error exceeds tolerance; increase mpol or "
            "crop inside the LCFS");
    for (int m = 0; m < input.mpol; ++m) {
        input.rbc_m.push_back(m);
        input.rbc_n.push_back(0);
        input.rbc_value.push_back(rc[m]);
        input.zbs_m.push_back(m);
        input.zbs_n.push_back(0);
        input.zbs_value.push_back(zs[m]);
        if (input.lasym) {
            input.rbs_m.push_back(m);
            input.rbs_n.push_back(0);
            input.rbs_value.push_back(rs[m]);
            input.zbc_m.push_back(m);
            input.zbc_n.push_back(0);
            input.zbc_value.push_back(zc[m]);
        }
    }
    input.rbcc = std::move(rc);
    input.zbsc = std::move(zs);
    input.rbss.assign(input.mpol, 0.0);
    input.zbcs.assign(input.mpol, 0.0);
    input.raxis_c = {result.lcfs.magnetic_axis.r};
    input.zaxis_s = {0.0};
    if (input.lasym) {
        input.rbsc = std::move(rs);
        input.zbcc = std::move(zc);
        input.rbcs.assign(input.mpol, 0.0);
        input.zbss.assign(input.mpol, 0.0);
        input.raxis_s = {0.0};
        input.zaxis_c = {result.lcfs.magnetic_axis.z};
    }
}
}  // namespace

LcfsShape calculate_lcfs_shape(std::span<const LcfsPoint> points,
                               double z_reference) {
    const auto polygon = clean_polygon(points);
    require(std::isfinite(z_reference), "LCFS Z reference must be finite");
    double r_min = polygon.front().r, r_max = r_min;
    double z_min = polygon.front().z, z_max = z_min;
    for (auto p : polygon) {
        r_min = std::min(r_min, p.r);
        r_max = std::max(r_max, p.r);
        z_min = std::min(z_min, p.z);
        z_max = std::max(z_max, p.z);
    }
    require(
        r_max > r_min && z_min < z_reference && z_reference < z_max,
        "LCFS needs nonzero width and points above and below the Z reference");
    const auto extremum_r = [&](double z) {
        double lo = std::numeric_limits<double>::infinity(), hi = -lo;
        for (auto p : polygon)
            if (p.z == z) {
                lo = std::min(lo, p.r);
                hi = std::max(hi, p.r);
            }
        return std::midpoint(lo, hi);
    };
    LcfsShape result;
    result.r_center = std::midpoint(r_min, r_max);
    result.minor_radius = 0.5 * (r_max - r_min);
    result.z_reference = z_reference;
    result.z_upper = z_max;
    result.z_lower = z_min;
    result.r_upper = extremum_r(z_max);
    result.r_lower = extremum_r(z_min);
    result.elongation_upper = (z_max - z_reference) / result.minor_radius;
    result.elongation_lower = (z_reference - z_min) / result.minor_radius;
    result.triangularity_upper =
        (result.r_center - result.r_upper) / result.minor_radius;
    result.triangularity_lower =
        (result.r_center - result.r_lower) / result.minor_radius;
    for (double value :
         {result.r_center, result.minor_radius, result.elongation_upper,
          result.elongation_lower, result.triangularity_upper,
          result.triangularity_lower})
        require(std::isfinite(value), "LCFS shape calculation overflow");
    return result;
}

LcfsShape calculate_lcfs_shape(const GfileLcfs& lcfs) {
    return calculate_lcfs_shape(lcfs.points, lcfs.magnetic_axis.z);
}
GfileLcfs read_gfile_lcfs(std::istream& input, std::optional<int> cocos) {
    return lcfs_of(load_gfile(input, cocos));
}
GfileLcfs read_gfile_lcfs(const std::filesystem::path& path,
                          std::optional<int> cocos) {
    auto input = open_gfile(path);
    return read_gfile_lcfs(input, cocos);
}

GfileInput gfile_to_cumes_input(std::istream& stream,
                                const GfileInputSettings& settings) {
    validate_settings(settings);
    auto data = load_gfile(stream, settings.cocos);
    const auto factors = data.cocos.factors_to<double>(Cocos(5));
    const double delta =
        factors.flux * (data.flux_LCFS - data.flux_magnetic_axis);
    require(std::isfinite(delta) && delta != 0.0,
            "G-file poloidal flux span must be finite and nonzero");
    std::vector<double> q = data.safety_factor;
    for (std::size_t k = 0; k < q.size(); ++k) {
        q[k] *= factors.safety_factor;
        const double f = factors.toroidal * data.f_pol[k];
        require(f != 0.0 && q[k] != 0.0 &&
                    std::signbit(f) == std::signbit(-delta * q[k]),
                "G-file q, F, and outward poloidal flux have inconsistent "
                "COCOS signs");
        require(data.pressure[k] >= 0.0, "G-file pressure must be nonnegative");
    }
    const FluxProfile flux(std::move(q), settings.psi_ratio);
    GfileInput result;
    result.lcfs = lcfs_of(data);
    result.psi_ratio = settings.psi_ratio;
    auto& input = result.input;
    input.mpol = settings.mpol;
    input.ntor = 0;
    input.nfp = 1;
    input.ntheta =
        settings.ntheta == 0 ? 2 * settings.mpol + 6 : settings.ntheta;
    input.nzeta = 1;
    input.ncurr = 0;
    input.delt = 0.9;
    input.phiedge = -PERIOD * delta * flux.edge();
    require(std::isfinite(input.phiedge) && input.phiedge != 0.0,
            "G-file integrated toroidal flux overflows or vanishes");
    input.pres_scale = 1.0;
    input.spres_ped = 1.0;
    input.bloat = 1.0;
    input.tcon0 = 1.0;
    input.curtor = factors.toroidal * data.current;
    input.pmass_type = "list";
    input.piota_type = "list";
    input.aphi = {1.0};
    input.ac = {0.0};
    input.stages = {settings.stage};
    for (std::size_t k = 0; k < settings.profile_points; ++k) {
        const double s = static_cast<double>(k) /
                         static_cast<double>(settings.profile_points - 1);
        const double x = flux.invert(s);
        input.am.push_back(FluxProfile::interpolate(data.pressure, x));
        input.ai.push_back(1.0 / flux.q(x));
    }
    double integral_iota = 0.0;
    for (std::size_t k = 1; k < input.ai.size(); ++k)
        integral_iota += 0.5 * (input.ai[k - 1] + input.ai[k]) /
                         static_cast<double>(input.ai.size() - 1);
    const double expected = settings.psi_ratio * delta;
    const double actual = -input.phiedge / PERIOD * integral_iota;
    result.poloidal_flux_relative_error =
        std::abs((actual - expected) / expected);
    require(std::isfinite(result.poloidal_flux_relative_error) &&
                result.poloidal_flux_relative_error <=
                    settings.poloidal_flux_tolerance,
            "resampled G-file iota changes poloidal flux beyond tolerance; "
            "increase profile_points");
    const auto sampled = sample_boundary(data, settings);
    fit_boundary(input, sampled, result, settings);
    const auto valid = validate_input_record(input);
    if (!valid) throw std::invalid_argument(valid.error());
    return result;
}
GfileInput gfile_to_cumes_input(const std::filesystem::path& path,
                                const GfileInputSettings& settings) {
    auto input = open_gfile(path);
    return gfile_to_cumes_input(input, settings);
}

std::string gfile_input_to_json(const GfileInput& converted) {
    const auto& p = converted.input;
    const auto valid = validate_input_record(p);
    if (!valid) throw std::invalid_argument(valid.error());
    require(p.ntor == 0 && p.nfp == 1 && p.ncurr == 0 && !p.lfreeb &&
                p.pmass_type == "list" && p.piota_type == "list" &&
                p.pcurr_type == "power_series" &&
                p.pflux_type == "power_series",
            "G-file JSON requires the generated fixed-boundary input model");
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::setprecision(17);
    out << "{\n  \"mpol\":" << p.mpol
        << ",\n  \"ntor\":0,\n  \"nfp\":1,\n  \"ntheta\":" << p.ntheta
        << ",\n  \"nzeta\":1,\n  \"ncurr\":0,\n  \"lfreeb\":false,\n  "
           "\"lasym\":"
        << (p.lasym ? "true" : "false") << ",\n  \"phiedge\":" << p.phiedge
        << ",\n  \"curtor\":" << p.curtor << ",\n  \"delt\":" << p.delt
        << ",\n  \"pres_scale\":" << p.pres_scale
        << ",\n  \"gamma\":" << p.adiabatic_index
        << ",\n  \"spres_ped\":" << p.spres_ped << ",\n  \"bloat\":" << p.bloat
        << ",\n  \"tcon0\":" << p.tcon0
        << ",\n  \"pmass_type\":\"list\",\n  \"piota_type\":\"list\"";
    const auto array = [&](const char* name, const auto& values) {
        out << ",\n  \"" << name << "\":[";
        for (std::size_t k = 0; k < values.size(); ++k) {
            if (k) out << ',';
            out << values[k];
        }
        out << ']';
    };
    array("am", p.am);
    array("ai", p.ai);
    array("ac", p.ac);
    array("aphi", p.aphi);
    array("raxis_c", p.raxis_c);
    array("zaxis_s", p.zaxis_s);
    if (p.lasym) {
        array("raxis_s", p.raxis_s);
        array("zaxis_c", p.zaxis_c);
    }
    std::vector<int> ns, iterations;
    std::vector<double> tolerance;
    for (const auto& stage : p.stages) {
        ns.push_back(stage.ns);
        iterations.push_back(stage.max_iter);
        tolerance.push_back(stage.ftol);
    }
    array("ns_array", ns);
    array("niter_array", iterations);
    array("ftol_array", tolerance);
    const auto boundary = [&](const char* name, const auto& m,
                              const auto& values) {
        out << ",\n  \"" << name << "\":[";
        for (std::size_t k = 0; k < m.size(); ++k) {
            if (k) out << ',';
            out << "{\"m\":" << m[k] << ",\"n\":0,\"value\":" << values[k]
                << '}';
        }
        out << ']';
    };
    boundary("rbc", p.rbc_m, p.rbc_value);
    boundary("zbs", p.zbs_m, p.zbs_value);
    if (p.lasym) {
        boundary("rbs", p.rbs_m, p.rbs_value);
        boundary("zbc", p.zbc_m, p.zbc_value);
    }
    out << "\n}\n";
    return out.str();
}

Status write_gfile_input(const std::string& path, const GfileInput& converted) {
    if (path.empty() || path.find('\0') != std::string::npos)
        return std::string("invalid G-file JSON output path");
    std::string temporary;
    bool owns_temporary = false;
    try {
        const auto json = gfile_input_to_json(converted);
        temporary = io_detail::temp_path_for(path);
        std::unique_ptr<std::FILE, CloseFile> file(
            std::fopen(temporary.c_str(), "wbx"));
        if (!file) throw std::runtime_error("cannot open G-file JSON output");
        owns_temporary = true;
        if (std::fwrite(json.data(), 1, json.size(), file.get()) != json.size())
            throw std::runtime_error("cannot write G-file JSON output");
        const auto error =
            io_detail::publish_atomic(file.release(), temporary, path);
        if (!error.empty()) return error;
        return {};
    } catch (const std::exception& error) {
        if (owns_temporary) std::remove(temporary.c_str());
        return std::string(error.what());
    }
}
}  // namespace equilibrium
