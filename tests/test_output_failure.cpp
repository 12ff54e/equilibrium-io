// Atomic publication and writer failure injection for every enabled codec.
#include "equilibrium/io.hpp"
#include "equilibrium/output_spec.hpp"
#include "equilibrium/run_report.hpp"
#include "equilibrium/writer.hpp"
#include "equilibrium/writer_helpers.hpp"
#include "hdf5_failure.hpp"
#include "input_fixture.hpp"
#include "test_support.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
using namespace equilibrium::test;

static bool file_exists(const char* path) {
    struct stat st;
    return stat(path, &st) == 0;
}

static void write_garbage(const char* path, const char* bytes, size_t n) {
    FILE* fp = fopen(path, "wb");
    if (!fp) {
        std::cerr << format("cannot seed {}\n", path);
        exit(1);
    }
    fwrite(bytes, 1, n, fp);
    fclose(fp);
}

template <typename T>
static bool write_via_host(const char* path, equilibrium::OutputFormat fmt) {
    equilibrium::EquilibriumSnapshot snapshot;
    snapshot.ns = 5;
    snapshot.mnmax = 2;
    for (auto& family : snapshot.families)
        family.assign(snapshot.family_size(), 0.0);
    equilibrium::RunReport report;
    report.build.scalar_type = sizeof(T) == sizeof(double) ? "double" : "float";
    report.input_params = input_fixture();
    equilibrium::StageReport stage;
    stage.ns = snapshot.ns;
    stage.effective_iterations = 3;
    stage.final_residual = {1e-10, 2e-10, 3e-10};
    report.stages.push_back(stage);
    auto writer = equilibrium::create_writer(fmt);
    return writer && writer
                         ->write_atomic(equilibrium::as_view(snapshot), report,
                                        {fmt, path})
                         .has_value();
}

template <typename T>
static void run_all() {
    std::cout << format("== {} precision ==\n",
                        sizeof(T) == sizeof(double) ? "double" : "float");

    // ---- open failure: path in a nonexistent directory ----
    const char* no_dir = "no_such_dir_cumes/state.bin";
    {
        // The binary writer must fail at fopen before reading the state.
        bool ok = write_via_host<T>(no_dir, equilibrium::OutputFormat::BINARY);
        check(!ok, "open failure: binary returns false");
        check(!file_exists(no_dir), "open failure: no partial file created");
    }
#ifdef EQUILIBRIUM_HAVE_NETCDF
    {
        bool ok = write_via_host<T>("no_such_dir_cumes/state.nc",
                                    equilibrium::OutputFormat::NETCDF);
        check(!ok, "open failure: netcdf returns false");
        check(!file_exists("no_such_dir_cumes/state.nc"),
              "open failure: netcdf no partial file");
    }
#endif
#ifdef EQUILIBRIUM_HAVE_HDF5
    {
        bool ok = write_via_host<T>("no_such_dir_cumes/state.h5",
                                    equilibrium::OutputFormat::HDF5);
        check(!ok, "open failure: hdf5 returns false");
        check(!file_exists("no_such_dir_cumes/state.h5"),
              "open failure: hdf5 no partial file");
    }
#endif

    // ---- atomic-rename failure: target is a directory with a valid suffix
    // ---- The dispatcher must route to the binary writer, so the target
    // carries a .bin suffix; the target IS a directory, so the writer writes
    // the temp fine but rename(temp, target) fails. It must return false,
    // remove the temp, and leave the directory untouched. Note: the temp path
    // is <target>.tmp.<pid>, which sits NEXT to the directory and writes fine.
    {
        const char* dir = "test_output_target_dir.bin";
        // Remove any prior run's leftover, then make a directory as the target.
        remove(dir);
        if (mkdir(dir, 0755) != 0) {
            std::cerr << format("cannot mkdir {}\n", dir);
            exit(1);
        }
        bool ok = write_via_host<T>(dir, equilibrium::OutputFormat::BINARY);
        check(!ok, "rename failure: target directory -> returns false");
        check(file_exists(dir), "rename failure: target directory untouched");
        // No stray temp file may remain next to the target.
        std::string tmp = std::string(dir) + ".tmp.";
        bool stray = false;
        // scan CWD for <dir>.tmp.* leftovers
        DIR* d = opendir(".");
        if (d) {
            struct dirent* e;
            while ((e = readdir(d))) {
                if (std::string(e->d_name).compare(0, tmp.size(), tmp) == 0)
                    stray = true;
            }
            closedir(d);
        }
        check(!stray, "rename failure: no stray temp file left");
        rmdir(dir);
    }

    // ---- truncation: a pre-existing file is clobbered, old bytes gone ----
    {
        const char* path = "test_output_trunc.bin";
        write_garbage(path, "GARBAGE-GARBAGE-GARBAGE-GARBAGE", 32);
        bool ok = write_via_host<T>(path, equilibrium::OutputFormat::BINARY);
        check(ok, "truncation: binary overwrite succeeds");
        // v1 header: magic (8) + version (4) + ns (4) + mnmax (4); the old
        // 'GARBAGE' bytes must be gone.
        FILE* fp = fopen(path, "rb");
        if (fp) {
            char magic[8] = {0};
            int version = 0, ns = 0, mnmax = 0;
            bool hdr = fread(magic, 1, 8, fp) == 8 &&
                       fread(&version, sizeof(int), 1, fp) == 1 &&
                       fread(&ns, sizeof(int), 1, fp) == 1 &&
                       fread(&mnmax, sizeof(int), 1, fp) == 1;
            fclose(fp);
            check(
                hdr && memcmp(magic, "CUMES001", 8) == 0 && version == 11 &&
                    ns == 5 && mnmax == 2,
                "truncation: binary header is the new run's (magic/ns/mnmax)");
        } else {
            check(false, "truncation: binary file readable");
        }
        remove(path);
    }
#ifdef EQUILIBRIUM_HAVE_NETCDF
    {
        const char* path = "test_output_trunc.nc";
        write_garbage(path, "GARBAGE-GARBAGE-GARBAGE-GARBAGE", 32);
        bool ok = write_via_host<T>(path, equilibrium::OutputFormat::NETCDF);
        check(ok, "truncation: netcdf overwrite succeeds");
        // NC_CLOBBER must have replaced the garbage; the file now starts with
        // the netCDF magic "CDF".
        FILE* fp = fopen(path, "rb");
        if (fp) {
            char magic[4] = {0};
            size_t got = fread(magic, 1, 3, fp);
            fclose(fp);
            check(got == 3 && memcmp(magic, "CDF", 3) == 0,
                  "truncation: netcdf file is a fresh CDF (not stale garbage)");
        } else {
            check(false, "truncation: netcdf file readable");
        }
        remove(path);
    }
#endif
#ifdef EQUILIBRIUM_HAVE_HDF5
    {
        const char* path = "test_output_trunc.h5";
        write_garbage(path, "GARBAGE-GARBAGE-GARBAGE-GARBAGE", 32);
        bool ok = write_via_host<T>(path, equilibrium::OutputFormat::HDF5);
        check(ok, "truncation: hdf5 overwrite succeeds");
        // HDF5 signature is "\211HDF\r\n\032\n".
        FILE* fp = fopen(path, "rb");
        if (fp) {
            const unsigned char sig[8] = {0x89, 'H',  'D',  'F',
                                          '\r', '\n', 0x1a, '\n'};
            unsigned char buf[8] = {0};
            size_t got = fread(buf, 1, 8, fp);
            fclose(fp);
            check(got == 8 && memcmp(buf, sig, 8) == 0,
                  "truncation: hdf5 file is a fresh HDF5 (not stale garbage)");
        } else {
            check(false, "truncation: hdf5 file readable");
        }
        remove(path);
    }
#endif
}

// ---------------------------------------------------------------------------
// Library-managed publication boundaries: the
// checked reopen/fsync/close/rename/directory-fsync chain publish_library_file
// implements for the NetCDF/HDF5 writers. Injected failures():
//
//   reopen  -> a missing temp fails before any rename;
//   fsync   -> /proc files reject fsync with EINVAL (a real, deterministic
//              fault on Linux; no fault-injection hook is needed);
//   rename  -> the target is a directory (covered end to end above too);
//   dir-fsync -> fsync_directory_of on /proc propagates the EINVAL instead of
//              ignoring it (the old helper ignored fsync AND close errors);
//
// close failures() cannot be injected portably without interposition; the
// checked close-after-fsync is the same three-line pattern as the fsync
// check and is exercised by the same code path here.
#if defined(__linux__) && !defined(__EMSCRIPTEN__)
static void run_publication_boundaries() {
    // Reopen boundary: the temp must exist and be readable.
    const char* dest1 = "test_output_pub_dest1.bin";
    remove(dest1);
    write_garbage(dest1, "OLD-DESTINATION", 15);
    {
        const std::string err = equilibrium::io_detail::publish_library_file(
            "test_output_pub_missing.tmp.x", dest1);
        check(!err.empty() && err.find("reopen") != std::string::npos,
              "library publish: missing temp fails at the reopen boundary");
        FILE* fp = fopen(dest1, "rb");
        char buf[16] = {0};
        const size_t got = fp ? fread(buf, 1, 15, fp) : 0;
        if (fp) fclose(fp);
        check(got == 15 && std::memcmp(buf, "OLD-DESTINATION", 15) == 0,
              "library publish: destination preserved on reopen failure");
    }
    remove(dest1);

    // fsync boundary: open succeeds, fsync fails (EINVAL on procfs).
    {
        const std::string err = equilibrium::io_detail::publish_library_file(
            "/proc/self/stat", "test_output_pub_dest2.bin");
        check(!err.empty() && err.find("fsync") != std::string::npos,
              "library publish: fsync failure propagates a typed reason");
        check(!file_exists("test_output_pub_dest2.bin"),
              "library publish: no destination after fsync failure");
    }

    // Directory-fsync boundary: checked and propagated (the old helper
    // ignored both the fsync and the close result).
    {
        const std::string err =
            equilibrium::io_detail::fsync_directory_of("/proc/self");
        check(!err.empty() && err.find("fsync") != std::string::npos,
              "directory fsync failure propagates a typed reason");
    }

    // Write/fflush boundary of the FILE* protocol: buffered writes to
    // /dev/full succeed until flush, where fflush fails — publish_atomic must
    // report it and remove nothing but its temp.
    {
        FILE* fp = fopen("/dev/full", "wb");
        check(fp != nullptr, "write fault: /dev/full opens");
        if (fp) {
            (void)fwrite("x", 1, 1, fp);  // buffered; fails at flush time
            const std::string err = equilibrium::io_detail::publish_atomic(
                fp, "/dev/full.tmp.x", "test_output_pub_dest3.bin");
            check(!err.empty() && err.find("fflush") != std::string::npos,
                  "write fault: fflush failure propagates a typed reason");
            check(!file_exists("test_output_pub_dest3.bin"),
                  "write fault: no destination after flush failure");
        }
    }
}

#endif

int main() {
    std::cout << "=== Output failure-injection matrix ===\n";
    run_all<double>();
    run_all<float>();
#if defined(__linux__) && !defined(__EMSCRIPTEN__)
    run_publication_boundaries();
#endif
#if defined(EQUILIBRIUM_HAVE_HDF5) && defined(__linux__) && \
    !defined(__EMSCRIPTEN__)
    try {
        hdf5_write_failures(
            "test_state_storage_failure.h5",
            [](const std::string& path) {
                require_hdf5(write_via_host<double>(
                                 path.c_str(), equilibrium::OutputFormat::HDF5),
                             "state HDF5 write failed");
            },
            [](const std::string& path) {
                // The reader probes optional datasets that may be absent.
                H5E_auto2_t handler = nullptr;
                void* handler_data = nullptr;
                require_hdf5(
                    H5Eget_auto2(H5E_DEFAULT, &handler, &handler_data) >= 0,
                    "cannot save HDF5 error handler");
                require_hdf5(H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr) >= 0,
                             "cannot silence optional HDF5 probes");
                const auto recovered = equilibrium::read(path);
                require_hdf5(
                    H5Eset_auto2(H5E_DEFAULT, handler, handler_data) >= 0,
                    "cannot restore HDF5 error handler");
                require_hdf5(recovered.has_value(),
                             "cannot read recovered state HDF5 file");
                const auto& state = recovered.value();
                require_hdf5(state.ns == 5 && state.mnmax == 2 &&
                                 state.families.size() == 6,
                             "recovered state dimensions changed");
                for (const auto& family : state.families)
                    require_hdf5(family == std::vector<double>(10, 0.0),
                                 "recovered state family changed");
            });
    } catch (const std::exception& error) { check(false, error.what()); }
#endif
    return summary();
}
