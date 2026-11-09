// test_checkpoint.cpp — versioned checkpoint.
//
// Verifies the Phase 2 "versioned checkpoint replacing CUMES_LOAD_INIT"
// deliverable: a versioned checkpoint round-trips and a corrupt
// magic/version/dimension is rejected.
#include "equilibrium/checkpoint.hpp"
#include "test_support.hpp"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>

#include <unistd.h>
using namespace equilibrium::test;

using equilibrium::EquilibriumSnapshot;

static std::string scratch(const char* name) {
    return std::string("test_checkpoint_") + name + "_" +
           std::to_string((long)getpid());
}

static EquilibriumSnapshot make_snapshot(int ns, int mnmax) {
    EquilibriumSnapshot s;
    s.ns = ns;
    s.mnmax = mnmax;
    const std::size_t n = s.family_size();
    for (int c = 0; c < 6; ++c) {
        s.families[c].resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            s.families[c][i] = c * 1000.0 + static_cast<double>(i) + 0.5;
        }
    }
    return s;
}

// A small but non-trivial embedded input record for the checkpoint trailer.
static equilibrium::InputParams make_params() {
    equilibrium::InputParams p;
    p.mpol = 2;
    p.ntor = 0;
    p.nfp = 1;
    p.ntheta = 10;
    p.nzeta = 1;
    p.ncurr = 0;
    p.delt = 0.9;
    p.phiedge = 1.0;
    p.am = {1.0};
    p.aphi = {1.0};
    p.raxis_c = p.zaxis_s = {0.0};
    p.rbc_m = {1, 0};
    p.rbc_n = {0, 0};
    p.rbc_value = {1.0, 2.0};
    p.zbs_m = {1};
    p.zbs_n = {0};
    p.zbs_value = {0.5};
    p.rbcc = {2.0, 1.0};
    p.rbss = {0.0, 0.0};
    p.zbsc = {0.0, 0.5};
    p.zbcs = {0.0, 0.0};
    p.stages.push_back({5, 100, 1e-12});
    return p;
}

static void test_checkpoint_round_trip() {
    EquilibriumSnapshot s = make_snapshot(5, 2);
    const equilibrium::InputParams params = make_params();
    const std::string path = scratch("rt");
    check(equilibrium::write_checkpoint(equilibrium::as_view(s), params, path)
              .has_value(),
          "checkpoint: write succeeds");
    equilibrium::InputParams back_params;
    auto back = equilibrium::read_binary(
        path, equilibrium::test::complete_request(true));
    if (back) back_params = *back.value().input;
    check(back.has_value() && equilibrium::test::same_state(s, back.value()),
          "checkpoint: round-trip preserves state");
    check(back_params == params,
          "checkpoint: round-trip preserves the input record");
    remove(path.c_str());
}

static void test_checkpoint_version_gate() {
    const std::string path = scratch("versions");
    for (int version : {1, 6, 7, 10}) {
        const auto snapshot = make_snapshot(4, 2);
        equilibrium::write_checkpoint(equilibrium::as_view(snapshot),
                                      make_params(), path);
        std::fstream file(path,
                          std::ios::in | std::ios::out | std::ios::binary);
        file.seekp(8);
        const std::int32_t value = version;
        file.write(reinterpret_cast<const char*>(&value), sizeof(value));
        file.close();
        check(!equilibrium::read_binary(path),
              "unsupported checkpoint version rejected");
    }
    remove(path.c_str());
}

static void test_checkpoint_rejection() {
    EquilibriumSnapshot s = make_snapshot(4, 2);
    const std::string path = scratch("bad");
    equilibrium::write_checkpoint(equilibrium::as_view(s), make_params(), path);
    // Corrupt the magic byte.
    {
        std::fstream f(path, std::ios::in | std::ios::out | std::ios::binary);
        f.seekp(0);
        f.put('X');
        f.close();
    }
    auto bad = equilibrium::read_binary(path);
    check(!bad.has_value(), "checkpoint: bad magic rejected");
    remove(path.c_str());

    // A truncated checkpoint (header only) is rejected.
    const std::string tpath = scratch("trunc");
    FILE* f = fopen(tpath.c_str(), "wb");
    fwrite("CUMECKP1", 1, 8, f);
    fclose(f);
    auto trunc = equilibrium::read_binary(tpath);
    check(!trunc.has_value(), "checkpoint: truncated header rejected");
    remove(tpath.c_str());
}

static void test_checkpoint_corrupt_huge_dimensions() {
    // A corrupt checkpoint header with ns*mnmax*48 bytes in [2^63, 2^64) used
    // to wrap the size_t -> long long cast in check_state_dimensions negative,
    // pass the file-size bound, and die in fam.resize(~2e17) with an uncaught
    // std::bad_alloc (std::terminate). The reader must return the
    // "dimensions implausible" error instead.
    const std::string path = scratch("hugedim");
    {
        FILE* f = fopen(path.c_str(), "wb");
        fwrite("CUMECKP1", 1, 8, f);
        const std::int32_t version = 9, precision = 0;
        const std::int32_t ns = 2147483647, mnmax = 95000000;
        fwrite(&version, sizeof(version), 1, f);
        fwrite(&precision, sizeof(precision), 1, f);
        fwrite(&ns, sizeof(ns), 1, f);
        fwrite(&mnmax, sizeof(mnmax), 1, f);
        const std::int32_t components = 6;
        fwrite(&components, sizeof(components), 1, f);
        fclose(f);
    }
    auto got = equilibrium::read_binary(path);
    check(!got.has_value(),
          "checkpoint: huge ns*mnmax rejected as implausible (no "
          "bad_alloc/terminate)");
    remove(path.c_str());
}

int main() {
    test_checkpoint_round_trip();
    test_checkpoint_version_gate();
    test_checkpoint_rejection();
    test_checkpoint_corrupt_huge_dimensions();
    return summary();
}
