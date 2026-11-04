#pragma once
#include "equilibrium/input_params.hpp"

namespace equilibrium::test {
inline InputParams input_fixture() {
    InputParams p;
    p.mpol = 2;
    p.ntor = 0;
    p.nfp = 1;
    p.ntheta = 10;
    p.nzeta = 1;
    p.delt = 0.9;
    p.phiedge = 1.0;
    p.am = {1.0};
    p.aphi = {1.0};
    p.raxis_c = {0.0};
    p.zaxis_s = {0.0};
    p.rbc_m = {1};
    p.rbc_n = {0};
    p.rbc_value = {1.0};
    p.zbs_m = {1};
    p.zbs_n = {0};
    p.zbs_value = {0.5};
    p.rbcc = {0.0, 1.0};
    p.rbss = {0.0, 0.0};
    p.zbsc = {0.0, 0.5};
    p.zbcs = {0.0, 0.0};
    p.stages = {{5, 100, 1e-12}};
    return p;
}
}  // namespace equilibrium::test
