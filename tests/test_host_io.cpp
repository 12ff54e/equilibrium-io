// test_host_io.cpp — host I/O: output-spec dispatch, versioned v1 round-trip
// and corruption rejection, and the writer failure matrix.
//
// This is the Phase 2 gate for the I/O half: the versioned v1 container
// round-trips its state (magic/version/header + six mode-major double families
// + provenance trailer) and rejects corruption, and every writer publishes
// atomically and fails cleanly on open/rename errors.
#include "equilibrium/checkpoint.hpp"
#include "equilibrium/output_spec.hpp"
#include "equilibrium/reader.hpp"
#include "equilibrium/run_report.hpp"
#include "equilibrium/writer.hpp"
#include "input_fixture.hpp"
#include "test_support.hpp"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>

#include <sys/stat.h>
#include <unistd.h>
using namespace equilibrium::test;

using equilibrium::EquilibriumSnapshot;
using equilibrium::OutputFormat;
using equilibrium::OutputSpec;
using equilibrium::RunReport;
using equilibrium::RunStatus;

static std::string scratch(const char* name) {
    return std::string("test_host_io_") + name + "_" +
           std::to_string((long)getpid());
}

// Deterministic snapshot: family c, index i -> value = c*1000 + i + 0.25.
static EquilibriumSnapshot make_snapshot(int ns, int mnmax) {
    EquilibriumSnapshot s;
    s.ns = ns;
    s.mnmax = mnmax;
    const std::size_t n = s.family_size();
    for (int c = 0; c < 6; ++c) {
        s.families[c].resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            s.families[c][i] = c * 1000.0 + static_cast<double>(i) + 0.25;
        }
    }
    s.ntheta = 3;
    s.nzeta = 2;
    for (int c = 0; c < EquilibriumSnapshot::HALF_FIELD_COUNT; ++c) {
        s.half_fields[c].resize(s.half_field_size());
        for (std::size_t i = 0; i < s.half_fields[c].size(); ++i) {
            s.half_fields[c][i] =
                10000.0 + c * 1000.0 + static_cast<double>(i) + 0.125;
        }
    }
    for (int c = 0; c < EquilibriumSnapshot::FULL_FIELD_COUNT; ++c) {
        s.full_fields[c].resize(s.full_field_size());
        for (std::size_t i = 0; i < s.full_fields[c].size(); ++i) {
            s.full_fields[c][i] =
                20000.0 + c * 1000.0 + static_cast<double>(i) + 0.375;
        }
    }
    return s;
}

static RunReport make_report() {
    RunReport r;
    r.status = RunStatus::CONVERGED;
    r.total_effective_iterations = 7;
    equilibrium::StageReport stage;
    stage.ns = 5;
    stage.effective_iterations = 7;
    stage.converged = true;
    stage.final_residual = {1e-14, 2e-15, 3e-16};
    r.stages.push_back(stage);
    r.build.revision = "deadbeef";
    r.build.build_type = "Release";
    r.build.scalar_type = "double";
    r.build.precision_policy = "verify-double";
    r.build.compile_flags = "";
    r.input.source_path = "inputs/solovev.json";
    r.runtime.gpu_name = "test-gpu";
    return r;
}

static void test_output_spec() {
    check(equilibrium::output_format_available(OutputFormat::BINARY),
          "output spec: binary always available");
    auto r = equilibrium::resolve_output_spec("state.bin");
    check(r.has_value() && r.value().format == OutputFormat::BINARY,
          "output spec: .bin resolves to binary");
    r = equilibrium::resolve_output_spec("state.H5");
    check(r.has_value() && r.value().format == OutputFormat::HDF5,
          "output spec: .H5 (case-insensitive) resolves to hdf5");
    r = equilibrium::resolve_output_spec("state.unknown");
    check(!r.has_value(), "output spec: unknown suffix rejected");
}

static void test_v1_round_trip() {
    EquilibriumSnapshot s = make_snapshot(4, 2);
    OutputSpec spec;
    spec.format = OutputFormat::BINARY;
    spec.path = scratch("v1").c_str();
    auto w = equilibrium::make_binary_writer();
    auto r = equilibrium::make_binary_reader();
    check(w != nullptr && r != nullptr, "v1: writer+reader factories");
    if (!w || !r) return;
    check(w->write_atomic(equilibrium::as_view(s), make_report(), spec)
              .has_value(),
          "v1: write succeeds");
    auto back = r->read(spec.path, equilibrium::test::state_request());
    check(back.has_value() && equilibrium::test::same_state(s, back.value()),
          "v1: round-trip preserves state");

    // Corrupt the magic -> reader rejects.
    {
        std::fstream f(spec.path,
                       std::ios::in | std::ios::out | std::ios::binary);
        f.seekp(0);
        f.put('X');
        f.close();
    }
    auto bad = r->read(spec.path, equilibrium::test::state_request());
    check(!bad.has_value(), "v1: bad magic rejected");
    remove(spec.path.c_str());
}

static void test_list_profiles_round_trip(bool lasym) {
    auto input = input_fixture();
    input.lasym = lasym;
    if (lasym) {
        input.raxis_s = {0.0};
        input.zaxis_c = {0.0};
        input.rbsc.assign(2, 0.0);
        input.rbcs.assign(2, 0.0);
        input.zbcc.assign(2, 0.0);
        input.zbss.assign(2, 0.0);
    }
    input.pmass_type = input.piota_type = input.pcurr_type = input.pflux_type =
        "list";
    input.am = {8.0, 2.0, 6.0};
    input.ai = {-0.2, 0.4};
    input.ac = {-1, 2};
    input.aphi = {1, 2};
    auto snapshot = make_snapshot(5, 2);
    if (lasym) snapshot.families.resize(12, snapshot.families.front());
    auto report = make_report();
    report.input_params = input;
    const std::string binary = scratch("list.bin");
    check(equilibrium::make_binary_writer()
              ->write_atomic(equilibrium::as_view(snapshot), report,
                             {OutputFormat::BINARY, binary})
              .has_value(),
          "list output: binary write");
    const auto decoded = equilibrium::make_binary_reader()->read(
        binary, equilibrium::test::complete_request());
    check(decoded && decoded.value().input == report.input_params,
          "list output: binary preserves complete input");
    remove(binary.c_str());
    const std::string checkpoint = scratch("list.checkpoint");
    check(equilibrium::write_checkpoint(as_view(snapshot), report.input_params,
                                        checkpoint)
              .has_value(),
          "list output: checkpoint write");
    const auto restored = equilibrium::read_binary(
        checkpoint, equilibrium::test::complete_request(true));
    check(restored && restored.value().input == report.input_params,
          "list output: checkpoint preserves complete input");
    remove(checkpoint.c_str());
}

static void test_corrupt_header_huge_dimensions() {
    // A corrupt v1 header with ns*mnmax*48 bytes in [2^63, 2^64) used to wrap
    // the size_t -> long long cast in check_state_dimensions negative, pass the
    // file-size bound, and die in fam.resize(~2e17) with an uncaught
    // std::bad_alloc (std::terminate). The reader must return the
    // "dimensions implausible" error instead.
    const std::string path = scratch("hugedim");
    {
        FILE* f = fopen(path.c_str(), "wb");
        const char magic[8] = {'C', 'U', 'M', 'E', 'S', '0', '0', '1'};
        fwrite(magic, 1, 8, f);
        const std::int32_t version = 1, ns = 2147483647, mnmax = 95000000;
        fwrite(&version, sizeof(version), 1, f);
        fwrite(&ns, sizeof(ns), 1, f);
        fwrite(&mnmax, sizeof(mnmax), 1, f);
        fclose(f);
    }
    auto r = equilibrium::make_binary_reader();
    auto got = r->read(path);
    check(!got.has_value(),
          "corrupt header: huge ns*mnmax rejected as implausible (no "
          "bad_alloc/terminate)");
    remove(path.c_str());
}

static void test_short_family_rejected() {
    // A snapshot with a family not sized ns*mnmax must fail the write BEFORE
    // touching the short vector's memory. The old pattern (`ok = false;` then
    // `ok = ok && write_f64_array(...)`) only skipped the OOB read through the
    // && short-circuit; the shared write_state_families now returns early on a
    // size mismatch, so no reordering can reintroduce the read.
    EquilibriumSnapshot s = make_snapshot(4, 2);
    // rmnss four elements short (32 bytes); shrink_to_fit keeps the capacity
    // exactly at the short size so any over-read crosses the allocation and is
    // visible to ASan/valgrind.
    s.families[3].resize(s.family_size() - 4);
    s.families[3].shrink_to_fit();
    OutputSpec spec;
    spec.format = OutputFormat::BINARY;
    spec.path = scratch("short");
    auto w = equilibrium::make_binary_writer();
    check(w->write_atomic(equilibrium::as_view(s), make_report(), spec)
                  .has_value() == false,
          "short family: write fails cleanly (no OOB read)");
    FILE* f = fopen(spec.path.c_str(), "rb");
    check(f == nullptr, "short family: no file published");
    if (f) fclose(f);
}

static void test_v1_unknown_precision_rejected() {
    // An unknown build scalar type must fail the v1 write (typed precision
    // tag), not silently record 0=double as the old string compare did.
    EquilibriumSnapshot s = make_snapshot(4, 2);
    OutputSpec spec;
    spec.format = OutputFormat::BINARY;
    spec.path = scratch("v1badprec");
    auto w = equilibrium::make_binary_writer();
    RunReport report = make_report();
    report.build.scalar_type = "single";
    check(!w->write_atomic(equilibrium::as_view(s), report, spec).has_value(),
          "v1: unknown precision tag rejected");
    FILE* f = fopen(spec.path.c_str(), "rb");
    check(f == nullptr, "v1: no file published on unknown precision tag");
    if (f) fclose(f);
}

static void test_failure_matrix() {
    EquilibriumSnapshot s = make_snapshot(4, 2);
    // open failure: a path in a nonexistent directory.
    OutputSpec spec;
    spec.format = OutputFormat::BINARY;
    spec.path = "no_such_dir/test_host_io_open.bin";
    auto w = equilibrium::make_binary_writer();
    check(w->write_atomic(equilibrium::as_view(s), make_report(), spec)
                  .has_value() == false,
          "failure: open failure returns false");

    // rename failure: target is a non-empty directory (temp writes fine,
    // rename() fails, temp removed, target untouched).
    const std::string dir = scratch("dir");
    mkdir(dir.c_str(), 0755);
    {
        FILE* f = fopen((dir + "/keep").c_str(), "w");
        if (f) {
            fputs("x", f);
            fclose(f);
        }
    }
    spec.path = dir;  // a directory, not a file
    check(w->write_atomic(equilibrium::as_view(s), make_report(), spec)
                  .has_value() == false,
          "failure: rename-over-directory returns false");
    // the directory and its contents are untouched
    FILE* keep = fopen((dir + "/keep").c_str(), "r");
    check(keep != nullptr, "failure: target directory untouched");
    if (keep) fclose(keep);
    remove((dir + "/keep").c_str());
    rmdir(dir.c_str());
}

int main() {
    test_output_spec();
    test_v1_round_trip();
    test_list_profiles_round_trip(false);
    test_list_profiles_round_trip(true);
    test_corrupt_header_huge_dimensions();
    test_short_family_rejected();
    test_v1_unknown_precision_rejected();
    test_failure_matrix();
    return summary();
}
