#include "equilibrium/io.hpp"

#include <stdexcept>
#include <string>

#include <emscripten.h>

namespace {

std::string error;

void convert(const std::string& input, const std::string& output) {
    const auto status = equilibrium::convert(input, output);
    if (!status) throw std::runtime_error(status.error());
}

}  // namespace

extern "C" EMSCRIPTEN_KEEPALIVE int equilibrium_format_convert(
    const char* input,
    const char* output) {
    try {
        error.clear();
        if (!input || !*input || !output || !*output)
            throw std::invalid_argument("missing conversion path");
        convert(input, output);
        return 0;
    } catch (const std::exception& caught) {
        error = caught.what();
        return 1;
    }
}

extern "C" EMSCRIPTEN_KEEPALIVE const char* equilibrium_format_error() {
    return error.c_str();
}
