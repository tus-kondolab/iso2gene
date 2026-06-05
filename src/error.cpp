#include "iso2gene/error.hpp"

namespace iso2gene {

Iso2GeneError::Iso2GeneError(ExitCode code, const std::string& message)
    : std::runtime_error(message), code_(code) {}

ExitCode Iso2GeneError::code() const noexcept {
    return code_;
}

} // namespace iso2gene
