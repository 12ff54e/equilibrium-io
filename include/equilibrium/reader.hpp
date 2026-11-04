// reader.hpp — host-memory reader interface.
#ifndef EQUILIBRIUM_INCLUDE_IO_READER_HPP_
#define EQUILIBRIUM_INCLUDE_IO_READER_HPP_

#include "equilibrium/equilibrium_snapshot.hpp"
#include "equilibrium/output_spec.hpp"
#include "equilibrium/result.hpp"
#include "equilibrium/run_report.hpp"

#include <functional>
#include <memory>
#include <optional>

namespace equilibrium {

class Reader {
   public:
    virtual ~Reader() = default;
    // Read spectral state and any scientific fields. When report is supplied,
    // decode and validate the provenance too. Binary state-only reads skip
    // the trailer; complete binary reads reject trailing data.
    virtual Result<EquilibriumSnapshot> read(
        const std::string& path,
        std::optional<std::reference_wrapper<RunReport>> report =
            std::nullopt) = 0;
};

// Full codec dispatch, linked through equilibrium::io. An unavailable codec
// returns nullptr. Use make_binary_reader with equilibrium::binary alone.
std::unique_ptr<Reader> make_reader(OutputFormat format);

// Binary-only factory (host I/O library; no backend defines required).
std::unique_ptr<Reader> make_binary_reader();

}  // namespace equilibrium

#endif  // EQUILIBRIUM_INCLUDE_IO_READER_HPP_
