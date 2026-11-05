#include <cstdio>
#include <filesystem>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <emscripten.h>

namespace {

using File = std::unique_ptr<std::FILE, decltype(&std::fclose)>;
std::vector<unsigned char> output;
std::string error;
int output_size = 0;

File open_file(const char* path, const char* mode) {
    if (!path || !*path) throw std::invalid_argument("missing file path");
    File file(std::fopen(path, mode), std::fclose);
    if (!file) throw std::runtime_error("could not open " + std::string(path));
    return file;
}

}  // namespace

extern "C" EMSCRIPTEN_KEEPALIVE int
equilibrium_file_write(const char* path, const unsigned char* bytes, int size) {
    try {
        error.clear();
        if (size < 0 || (size && !bytes))
            throw std::invalid_argument("invalid file data");
        if (!path || !*path) throw std::invalid_argument("missing file path");
        const auto parent = std::filesystem::path(path).parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);
        auto file = open_file(path, "wb");
        if (std::fwrite(bytes, 1, static_cast<std::size_t>(size), file.get()) !=
                static_cast<std::size_t>(size) ||
            std::fflush(file.get()))
            throw std::runtime_error("could not write " + std::string(path));
        return 0;
    } catch (const std::exception& caught) {
        error = caught.what();
        return 1;
    }
}

extern "C" EMSCRIPTEN_KEEPALIVE const unsigned char* equilibrium_file_read(
    const char* path) {
    try {
        error.clear();
        output_size = 0;
        auto file = open_file(path, "rb");
        if (std::fseek(file.get(), 0, SEEK_END))
            throw std::runtime_error("could not seek input");
        const auto size = std::ftell(file.get());
        if (size < 0 || size > std::numeric_limits<int>::max())
            throw std::runtime_error("file exceeds Wasm input limit");
        std::rewind(file.get());
        output.resize(static_cast<std::size_t>(size) + 1);
        if (std::fread(output.data(), 1, static_cast<std::size_t>(size),
                       file.get()) != static_cast<std::size_t>(size))
            throw std::runtime_error("could not read " + std::string(path));
        output_size = static_cast<int>(size);
        return output.data();
    } catch (const std::exception& caught) {
        error = caught.what();
        return nullptr;
    }
}

extern "C" EMSCRIPTEN_KEEPALIVE int equilibrium_file_size() {
    return output_size;
}

extern "C" EMSCRIPTEN_KEEPALIVE int equilibrium_file_remove(const char* path) {
    return path ? std::remove(path) : -1;
}

extern "C" EMSCRIPTEN_KEEPALIVE const char* equilibrium_file_error() {
    return error.c_str();
}
