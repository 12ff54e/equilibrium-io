// reader.hpp — host-memory reader interface.
#ifndef EQUILIBRIUM_INCLUDE_IO_READER_HPP_
#define EQUILIBRIUM_INCLUDE_IO_READER_HPP_

#include "equilibrium/equilibrium.hpp"
#include "equilibrium/output_spec.hpp"
#include "equilibrium/result.hpp"

#include <memory>

namespace equilibrium {

class Reader {
   public:
    virtual ~Reader() = default;
    virtual Result<Equilibrium> read(
        const std::string& path,
        const ReadRequest& request = ReadRequest::all()) = 0;
};

// Full codec dispatch, linked through equilibrium::io. An unavailable codec
// returns nullptr. Use make_binary_reader with equilibrium::binary alone.
std::unique_ptr<Reader> make_reader(OutputFormat format);

// Binary-only factory (host I/O library; no backend defines required).
std::unique_ptr<Reader> make_binary_reader();

}  // namespace equilibrium

#endif  // EQUILIBRIUM_INCLUDE_IO_READER_HPP_
