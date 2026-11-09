// Complete latest-schema round trips in every enabled codec.
#include "equilibrium/checkpoint.hpp"
#include "equilibrium/output_spec.hpp"
#include "equilibrium/reader.hpp"
#include "equilibrium/writer.hpp"
#include "input_fixture.hpp"
#include "test_support.hpp"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include <unistd.h>
using namespace equilibrium::test;

using equilibrium::EquilibriumSnapshot;
using equilibrium::OutputFormat;
using equilibrium::OutputSpec;
using equilibrium::RunReport;

static std::string scratch(const char* name) {
    return std::string("test_io_golden_") + name + "_" +
           std::to_string((long)getpid());
}

// Deterministic pattern: family c, index i -> value = c*1000 + i + 0.25
// (0.25 is exactly representable in both float and double, so the T->double
// conversion is lossless in both precisions).
template <typename T>
static EquilibriumSnapshot make_snapshot(int ns, int mnmax) {
    EquilibriumSnapshot snapshot;
    snapshot.ns = ns;
    snapshot.mnmax = mnmax;
    for (int c = 0; c < 6; ++c) {
        auto& family = snapshot.families[c];
        family.resize(snapshot.family_size());
        for (std::size_t i = 0; i < family.size(); ++i)
            family[i] =
                static_cast<T>(c * 1000.0 + static_cast<double>(i) + 0.25);
    }
    return snapshot;
}

static void add_derived_fields(EquilibriumSnapshot& snapshot) {
    snapshot.ntheta = 3;
    snapshot.nzeta = 2;
    for (int c = 0; c < EquilibriumSnapshot::HALF_FIELD_COUNT; ++c) {
        snapshot.half_fields[c].resize(snapshot.half_field_size());
        for (std::size_t i = 0; i < snapshot.half_fields[c].size(); ++i) {
            snapshot.half_fields[c][i] =
                10000.0 + c * 1000.0 + static_cast<double>(i) + 0.125;
        }
    }
    for (int c = 0; c < EquilibriumSnapshot::FULL_FIELD_COUNT; ++c) {
        snapshot.full_fields[c].resize(snapshot.full_field_size());
        for (std::size_t i = 0; i < snapshot.full_fields[c].size(); ++i) {
            snapshot.full_fields[c][i] =
                20000.0 + c * 1000.0 + static_cast<double>(i) + 0.375;
        }
    }
}

// Complete provenance and saved input must survive each container exactly.
static RunReport make_io_report(const equilibrium::InputParams& input) {
    RunReport r;
    r.status = equilibrium::RunStatus::CONVERGED;
    r.total_effective_iterations = 21;
    equilibrium::StageReport s1;
    s1.ns = 5;
    s1.effective_iterations = 9;
    s1.converged = true;
    s1.final_residual = {1e-10, 2e-11, 3e-12};
    s1.restarts = {equilibrium::RestartEvent{2}, equilibrium::RestartEvent{5}};
    equilibrium::StageReport s2;
    s2.ns = 8;
    s2.effective_iterations = 12;
    s2.converged = true;
    s2.final_residual = {4e-13, 5e-14, 6e-15};
    s2.restarts = {equilibrium::RestartEvent{3}};
    r.stages = {s1, s2};
    r.build.revision = "abc123";
    r.build.dirty = true;
    r.build.build_type = "Release";
    r.build.scalar_type = "double";
    r.build.precision_policy = "verify-double";
    r.build.compile_flags = "";
    r.input.source_path = "inputs/w7x.json";
    r.input.source_hash = "0123456789abcdef";
    r.runtime.gpu_name = "NVIDIA TITAN Xp";
    r.runtime.driver = "580.173.02";
    r.runtime.runtime = "12010";
    r.runtime.toolkit = "12.1";
    r.input_params = input;
    r.input_params.embedded_makegrid_parameters =
        equilibrium::MakegridParametersSpec{false, false, 1,   2.0, 6.0,
                                            201,   -2.0,  2.0, 201, 1};
    return r;
}

static bool reports_equal(const RunReport& a, const RunReport& b) {
    if (a.status != b.status ||
        a.total_effective_iterations != b.total_effective_iterations) {
        return false;
    }
    if (a.build.revision != b.build.revision ||
        a.build.dirty != b.build.dirty ||
        a.build.build_type != b.build.build_type ||
        a.build.scalar_type != b.build.scalar_type ||
        a.build.precision_policy != b.build.precision_policy ||
        a.build.compile_flags != b.build.compile_flags ||
        a.input.source_path != b.input.source_path ||
        a.input.source_hash != b.input.source_hash ||
        a.runtime.gpu_name != b.runtime.gpu_name ||
        a.runtime.driver != b.runtime.driver ||
        a.runtime.runtime != b.runtime.runtime ||
        a.runtime.toolkit != b.runtime.toolkit) {
        return false;
    }
    if (!(a.input_params == b.input_params)) return false;
    if (a.stages.size() != b.stages.size()) return false;
    for (size_t g = 0; g < a.stages.size(); ++g) {
        const auto& x = a.stages[g];
        const auto& y = b.stages[g];
        if (x.ns != y.ns || x.effective_iterations != y.effective_iterations ||
            x.converged != y.converged ||
            x.final_residual.fsqr != y.final_residual.fsqr ||
            x.final_residual.fsqz != y.final_residual.fsqz ||
            x.final_residual.fsql != y.final_residual.fsql ||
            x.restarts.size() != y.restarts.size()) {
            return false;
        }
        for (size_t k = 0; k < x.restarts.size(); ++k) {
            if (x.restarts[k].iteration != y.restarts[k].iteration)
                return false;
        }
    }
    return true;
}

// State and complete run metadata must survive every codec.
template <typename T>
static void check_latest_round_trip(const EquilibriumSnapshot& snap,
                                    const equilibrium::InputParams& input,
                                    OutputFormat fmt,
                                    const char* tag) {
    const std::string path = scratch(tag);
    OutputSpec spec;
    spec.format = fmt;
    spec.path = path;
    auto w = equilibrium::create_writer(spec.format);
    auto r = equilibrium::make_reader(spec.format);
    check(w != nullptr && r != nullptr, "latest factories");
    if (w && r) {
        RunReport report = make_io_report(input);
        report.build.scalar_type =
            sizeof(T) == sizeof(double) ? "double" : "float";
        check(w->write_atomic(equilibrium::as_view(snap), report, spec)
                  .has_value(),
              "latest write succeeds");
        RunReport back;
        auto snap_back = r->read(path, equilibrium::test::complete_request());
        if (snap_back) back = equilibrium::test::report_from(snap_back.value());
        check(snap_back.has_value() &&
                  equilibrium::test::same_state(snap, snap_back.value()),
              "latest state round trip");
        check(reports_equal(report, back),
              "latest complete RunReport + restart metadata round trip");
    }
    remove(path.c_str());
}

template <typename T>
static void run_precision() {
    std::cout << format("== {} precision ==\n",
                        sizeof(T) == sizeof(double) ? "double" : "float");
    const int ns = 5, mnmax = 2;
    const auto state = make_snapshot<T>(ns, mnmax);

    auto input = input_fixture();
    input.ntheta = 3;
    input.nzeta = 2;
    // ---- latest writer round-trips the snapshot + provenance trailer --
    {
        const EquilibriumSnapshot snap = state;
        const std::string latestPath = scratch("latest");
        OutputSpec spec;
        spec.format = OutputFormat::BINARY;
        spec.path = latestPath;
        auto w = equilibrium::create_writer(spec.format);
        auto r = equilibrium::make_reader(spec.format);
        check(w != nullptr && r != nullptr, "latest: writer+reader factories");
        if (w && r) {
            // Read back WITH the report: the provenance trailer
            // (precision-policy fields + stage records) must round-trip, not
            // just the state. (The trailer sequences were once out of sync and
            // a state-only read could not see it — see docs/output-formats.md.)
            RunReport report = make_io_report(input);
            report.build.scalar_type =
                sizeof(T) == sizeof(double) ? "double" : "float";
            check(w->write_atomic(equilibrium::as_view(snap), report, spec)
                      .has_value(),
                  "latest: write succeeds");
            RunReport back;
            auto snap_back =
                r->read(latestPath, equilibrium::test::complete_request());
            if (snap_back)
                back = equilibrium::test::report_from(snap_back.value());
            check(snap_back.has_value() &&
                      equilibrium::test::same_state(snap, snap_back.value()),
                  "latest: round-trip preserves state");
            check(reports_equal(report, back),
                  "latest: round-trip preserves the provenance trailer");
        }
        remove(latestPath.c_str());
    }

    // ---- NetCDF/HDF5 containers ----------------------------------------
    // latest: the complete RunReport + restart metadata round trip.
    EquilibriumSnapshot snap2 = state;
    add_derived_fields(snap2);
#ifdef EQUILIBRIUM_HAVE_NETCDF
    check_latest_round_trip<T>(snap2, input, OutputFormat::NETCDF, "latestnc");
#endif
#ifdef EQUILIBRIUM_HAVE_HDF5
    check_latest_round_trip<T>(snap2, input, OutputFormat::HDF5, "latesth5");
#endif

    // ---- versioned checkpoint round-trips the snapshot ------------
    {
        const EquilibriumSnapshot snap = state;
        const std::string ckpt = scratch("ckpt");
        const equilibrium::InputParams ip = input;
        check(
            equilibrium::write_checkpoint(equilibrium::as_view(snap), ip, ckpt)
                .has_value(),
            "checkpoint: write succeeds");
        equilibrium::InputParams back_params;
        auto back = equilibrium::read_binary(
            ckpt, equilibrium::test::complete_request(true));
        if (back) back_params = *back.value().input;
        check(back.has_value() &&
                  equilibrium::test::same_state(snap, back.value()),
              "checkpoint: round-trip preserves state");
        check(back_params == ip,
              "checkpoint: round-trip preserves the input record");
        remove(ckpt.c_str());
    }
}

int main() {
    std::cout << "=== latest writer round-trip gate ===\n";
    run_precision<double>();
    run_precision<float>();
    return summary();
}
