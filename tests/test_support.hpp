#ifndef EQUILIBRIUM_TESTS_TEST_H_
#define EQUILIBRIUM_TESTS_TEST_H_
// Shared C++20 assertions and fixture directories.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <span>
#include <string_view>

#include <unistd.h>

// Formatting also supports GCC 12, which lacks std::format.
#include "test_format.hpp"

namespace equilibrium::test {

// Each CTest executable has its own working directory in the build tree.
class TempDir {
   public:
    TempDir() : dir_("fixtures_" + std::to_string(getpid())) {
        std::filesystem::create_directories(dir_);
    }
    ~TempDir() { std::filesystem::remove_all(dir_); }
    const std::string& path() const { return dir_; }
    bool ok() const { return std::filesystem::is_directory(dir_); }

   private:
    std::string dir_;
};

// The single shared failure counter (function-local static, so header-only and
// shared across TUs). Capture `int before = failures();` to count only the
// failures a sub-check adds.
inline int& failures() {
    static int n = 0;
    return n;
}

// Assert `ok`; on failure print "FAIL <msg>" and increment the counter.
inline void check(bool ok, std::string_view msg) {
    if (ok) {
        std::cout << "PASS " << msg << '\n';
    } else {
        std::cout << "FAIL " << msg << '\n';
        ++failures();
    }
}

// Print the final summary and return the process exit code (0 = all pass).
inline int summary() {
    if (failures()) {
        std::cout << failures() << " FAILURES\n";
        return 1;
    }
    std::cout << "ALL PASS\n";
    return 0;
}

}  // namespace equilibrium::test

#endif  // EQUILIBRIUM_TESTS_TEST_H_
