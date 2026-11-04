#include "equilibrium/io.hpp"

#include <iostream>

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: equilibrium-convert INPUT OUTPUT.{bin,nc,h5}\n";
        return 2;
    }
    const auto status = equilibrium::convert(argv[1], argv[2]);
    if (!status) {
        std::cerr << status.error() << '\n';
        return 1;
    }
}
