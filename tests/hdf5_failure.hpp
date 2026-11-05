#pragma once

#if defined(EQUILIBRIUM_HAVE_HDF5) && defined(__linux__) && \
    !defined(__EMSCRIPTEN__)

#include <cerrno>
#include <csignal>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <hdf5.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

namespace equilibrium::test {

inline void require_hdf5(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

inline std::size_t open_descriptors() {
    std::size_t count = 0;
    for ([[maybe_unused]] const auto& entry :
         std::filesystem::directory_iterator("/proc/self/fd"))
        ++count;
    return count;
}

// Run each persistent storage fault in a child, including HDF5 shutdown and
// normal exit. Drain stderr through a pipe: RLIMIT_FSIZE would truncate a log
// file too, and H5close can return success while printing a shutdown warning.
template <class Write, class Verify>
void hdf5_write_failures(const std::string& path, Write write, Verify verify) {
    const std::string sentinel = "previous destination";
    for (const rlim_t bytes : {0, 256, 4096}) {
        {
            std::ofstream seed(path, std::ios::binary);
            seed << sentinel;
            require_hdf5(seed.good(), "cannot seed HDF5 failure destination");
        }
        int pipe_fd[2];
        require_hdf5(pipe(pipe_fd) == 0, "cannot create HDF5 diagnostic pipe");
        std::cout.flush();
        std::cerr.flush();
        const pid_t child = fork();
        if (child < 0) {
            close(pipe_fd[0]);
            close(pipe_fd[1]);
            throw std::runtime_error("cannot fork HDF5 failure test");
        }
        if (child == 0) {
            close(pipe_fd[0]);
            if (dup2(pipe_fd[1], STDERR_FILENO) < 0) std::exit(2);
            close(pipe_fd[1]);
            try {
                rlimit saved{};
                require_hdf5(getrlimit(RLIMIT_FSIZE, &saved) == 0,
                             "cannot read file-size limit");
                require_hdf5(saved.rlim_max >= bytes,
                             "hard file-size limit is too small for test");
                const auto other_path = path + ".other.h5";
                write(other_path);
                const hid_t other =
                    H5Fopen(other_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
                require_hdf5(other >= 0, "cannot open unrelated HDF5 file");
                const auto objects = H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL);
                const auto descriptors = open_descriptors();
                require_hdf5(std::signal(SIGXFSZ, SIG_IGN) != SIG_ERR,
                             "cannot ignore file-size signal");
                const rlimit limited{bytes, saved.rlim_max};
                require_hdf5(setrlimit(RLIMIT_FSIZE, &limited) == 0,
                             "cannot limit file size");
                for (int repeat = 0; repeat < 2; ++repeat) {
                    bool failed = false;
                    try {
                        write(path);
                    } catch (const std::exception&) { failed = true; }
                    require_hdf5(failed, "HDF5 storage failure not reported");
                    require_hdf5(
                        H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL) == objects,
                        "HDF5 failure leaked object identifiers");
                    require_hdf5(open_descriptors() == descriptors,
                                 "HDF5 failure leaked file descriptors");
                    require_hdf5(H5Iis_valid(other) > 0,
                                 "HDF5 failure closed an unrelated file");
                    std::ifstream previous(path, std::ios::binary);
                    require_hdf5(
                        std::string(std::istreambuf_iterator<char>(previous),
                                    {}) == sentinel,
                        "HDF5 failure replaced destination");
                    for (const auto& entry :
                         std::filesystem::directory_iterator("."))
                        require_hdf5(
                            !entry.path().filename().string().starts_with(
                                path + ".tmp."),
                            "HDF5 failure left a temporary file");
                }
                require_hdf5(H5Fclose(other) >= 0,
                             "cannot close unrelated HDF5 file");
                std::filesystem::remove(other_path);
                // The storage fault is still active during shutdown.
                require_hdf5(H5close() >= 0, "HDF5 shutdown failed");
                require_hdf5(setrlimit(RLIMIT_FSIZE, &saved) == 0,
                             "cannot restore file-size limit");
                write(path);
                verify(path);
            } catch (const std::exception& error) {
                std::cerr << error.what() << '\n';
                std::exit(1);
            }
            std::exit(0);  // Exercise HDF5's atexit handler and sanitizers.
        }
        close(pipe_fd[1]);
        std::string diagnostics;
        char buffer[4096];
        ssize_t count;
        while ((count = ::read(pipe_fd[0], buffer, sizeof(buffer))) != 0) {
            if (count < 0) {
                if (errno == EINTR) continue;
                break;
            }
            diagnostics.append(buffer, static_cast<std::size_t>(count));
        }
        close(pipe_fd[0]);
        int status = 0;
        pid_t waited;
        do {
            waited = waitpid(child, &status, 0);
        } while (waited < 0 && errno == EINTR);
        if (count < 0 || waited != child || !WIFEXITED(status) ||
            WEXITSTATUS(status) != 0 || !diagnostics.empty())
            throw std::runtime_error("HDF5 failure/shutdown test at " +
                                     std::to_string(bytes) +
                                     " bytes: " + diagnostics);
        std::filesystem::remove(path);
    }
}

}  // namespace equilibrium::test

#endif
