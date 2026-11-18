#pragma once
#include <array>
#include <cmath>
#include <iomanip>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

namespace gfile_test {
constexpr double R0 = 3.0, A = 0.6, B = 0.84, PSI = 0.15;
constexpr double Q0 = 1.7, Q1 = 0.8, P0 = 18000.0;
struct Settings {
    int cocos = 5;
    int current_sign = 1;
    int field_sign = 1;
    double axis_z = 0.2;
    double tilt = 0.0;
    double grid_flux_scale = 1.0;
    bool zero_profiles = false;
    bool degenerate_angles = false;
    bool flip_q = false;
    bool reverse_boundary = false;
    bool close_boundary = false;
    bool nonstar = false;
    bool fortran_d = false;
    std::string header;
};
// Independent Table I signs, not Cocos::factors_to from the implementation.
inline constexpr std::array<std::array<int, 3>, 8> COCOS{{{1, 1, 1},
                                                          {1, -1, 1},
                                                          {-1, 1, -1},
                                                          {-1, -1, -1},
                                                          {1, 1, -1},
                                                          {1, -1, -1},
                                                          {-1, 1, 1},
                                                          {-1, -1, 1}}};
inline std::string fixture(const Settings& settings = {}) {
    const auto [bp, phi, theta] = COCOS[(settings.cocos % 10) - 1];
    const double flux_factor =
        phi * bp * (settings.cocos > 10 ? 2 * std::numbers::pi : 1);
    constexpr int N = 17, NB = 128;
    const int ip = settings.current_sign, bt = settings.field_sign;
    std::ostringstream out;
    out << std::left << std::setw(48)
        << (settings.header.empty()
                ? "analytic COCOS=" + std::to_string(settings.cocos)
                : settings.header)
        << std::right << std::setw(4) << 0 << std::setw(4) << N << std::setw(4)
        << N << '\n';
    const auto header_size = out.str().size();
    out << std::scientific << std::uppercase << std::setprecision(9);
    const auto values = [&](const std::vector<double>& data) {
        for (std::size_t i = 0; i < data.size(); ++i) {
            out << std::setw(16) << data[i];
            if (i % 5 == 4 || i + 1 == data.size()) out << '\n';
        }
    };
    values({2.4, 2.8, R0, 1.8, settings.axis_z});
    values({R0, settings.axis_z, flux_factor * ip * .7,
            flux_factor * ip * (.7 + PSI), phi * bt * 2.0});
    values({phi * ip * 2.e6, 0, 0, 0, 0});
    values({0, 0, 0, 0, 0});
    std::vector<double> f, p, q, pp(N, -P0 * .5 / (flux_factor * ip * PSI));
    for (int i = 0; i < N; ++i) {
        const double u = static_cast<double>(i) / (N - 1);
        f.push_back(settings.zero_profiles ? 0.0 : phi * bt * 6.0);
        p.push_back(P0 * (1 - .5 * u));
        q.push_back(settings.zero_profiles ? 0.0
                                           : theta * ip * bt * (Q0 + Q1 * u));
    }
    if (settings.flip_q) q[N / 2] = -q[N / 2];
    values(f);
    values(p);
    values(std::vector<double>(N));
    values(pp);
    std::vector<double> flux;
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) {
            const double r = -1.2 + 2.4 * i / (N - 1);
            const double z = -1.4 + 2.8 * j / (N - 1);
            flux.push_back(flux_factor * ip *
                           (.7 + settings.grid_flux_scale * PSI *
                                     (r * r / (A * A) + z * z / (B * B) +
                                      2 * settings.tilt * r * z / (A * B))));
        }
    values(flux);
    values(q);
    if (settings.degenerate_angles) {
        out << "3 0\n";
        values({4, -1e-17, 5, -6e-17, 6, -1.5e-16});
        return out.str();
    }
    out << (NB + (settings.close_boundary ? 1 : 0)) << " 0\n";
    std::vector<double> points;
    for (int k = 0; k < NB + (settings.close_boundary ? 1 : 0); ++k) {
        int index = settings.reverse_boundary ? NB - k : k;
        if (settings.nonstar && k == 2) index = 20;
        const double t = 2 * std::numbers::pi * index / NB;
        const double c = std::cos(t), s = std::sin(t);
        const double radius = A * B /
                              std::sqrt(B * B * c * c + A * A * s * s +
                                        2 * settings.tilt * A * B * c * s);
        points.push_back(R0 + radius * c);
        points.push_back(settings.axis_z + radius * s);
    }
    values(points);
    auto text = out.str();
    if (settings.fortran_d)
        for (std::size_t k = header_size; k < text.size(); ++k)
            if (text[k] == 'E') text[k] = 'D';
    return text;
}
inline double poloidal_fraction(double s, double crop = 1) {
    const double target = s * (Q0 * crop + .5 * Q1 * crop * crop);
    return 2 * target / (Q0 + std::sqrt(Q0 * Q0 + 2 * Q1 * target));
}
}  // namespace gfile_test
