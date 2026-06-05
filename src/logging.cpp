#include "iso2gene/logging.hpp"

#include <fstream>
#include <iostream>

#include "iso2gene/error.hpp"

namespace iso2gene {

void Logger::info(const std::string& message) {
    std::cerr << "INFO: " << message << '\n';
}

void Logger::warn(const std::string& message) {
    warnings_.push_back(message);
    std::cerr << "WARNING: " << message << '\n';
}

bool Logger::has_warnings() const {
    return !warnings_.empty();
}

const std::vector<std::string>& Logger::warnings() const {
    return warnings_;
}

void Logger::write_warnings(const std::string& path) const {
    std::ofstream out(path);
    if (!out) {
        throw Iso2GeneError(ExitCode::io_error, "failed to open warnings file for writing: " + path);
    }
    for (const std::string& warning : warnings_) {
        out << warning << '\n';
    }
}

} // namespace iso2gene
