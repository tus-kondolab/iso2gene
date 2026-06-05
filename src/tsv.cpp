#include "iso2gene/tsv.hpp"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include "iso2gene/error.hpp"

namespace iso2gene {

std::vector<std::string> split_tsv_line(const std::string& line) {
    std::string cleaned = line;
    // Defensive for direct callers; TextReader already normalizes CRLF/CR line endings.
    if (!cleaned.empty() && cleaned.back() == '\r') {
        cleaned.pop_back();
    }

    std::vector<std::string> fields;
    std::string::size_type start = 0;
    while (true) {
        const std::string::size_type tab = cleaned.find('\t', start);
        if (tab == std::string::npos) {
            fields.push_back(cleaned.substr(start));
            break;
        }
        fields.push_back(cleaned.substr(start, tab - start));
        start = tab + 1;
    }
    return fields;
}

double parse_double_strict(
    const std::string& value,
    const std::string& path,
    std::size_t line_number,
    const std::string& column_name
) {
    if (value.empty()) {
        std::ostringstream msg;
        msg << path << ":" << line_number << ": empty numeric field '" << column_name << "'";
        throw Iso2GeneError(ExitCode::parse_error, msg.str());
    }

    errno = 0;
    char* end = nullptr;
    const double parsed = std::strtod(value.c_str(), &end);
    if (end == value.c_str() || *end != '\0' || errno == ERANGE || !std::isfinite(parsed)) {
        std::ostringstream msg;
        msg << path << ":" << line_number << ": invalid numeric field '" << column_name
            << "': " << value;
        throw Iso2GeneError(ExitCode::parse_error, msg.str());
    }
    return parsed;
}

TsvReader::TsvReader(const std::string& path)
    : reader_(path) {}

bool TsvReader::read_row(std::vector<std::string>& fields) {
    std::string line;
    if (!reader_.read_line(line)) {
        return false;
    }
    fields = split_tsv_line(line);
    return true;
}

std::size_t TsvReader::line_number() const noexcept {
    return reader_.line_number();
}

const std::string& TsvReader::path() const noexcept {
    return reader_.path();
}

int find_column(
    const std::vector<std::string>& header,
    const std::string& column_name,
    const std::string& path
) {
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (header[i] == column_name) {
            return static_cast<int>(i);
        }
    }
    throw Iso2GeneError(
        ExitCode::parse_error,
        path + ": missing required column '" + column_name + "'"
    );
}

} // namespace iso2gene
