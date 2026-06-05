#pragma once

#include <stdexcept>
#include <string>

namespace iso2gene {

enum class ExitCode {
    success = 0,
    input_error = 1,
    io_error = 2,
    parse_error = 3,
    internal_error = 4
};

class Iso2GeneError : public std::runtime_error {
public:
    Iso2GeneError(ExitCode code, const std::string& message);

    ExitCode code() const noexcept;

private:
    ExitCode code_;
};

} // namespace iso2gene
