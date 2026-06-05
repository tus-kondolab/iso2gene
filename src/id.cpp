#include "iso2gene/id.hpp"

namespace iso2gene {

std::string normalize_transcript_id(const std::string& id, const IdOptions& options) {
    std::string normalized = id;

    if (options.ignore_after_bar) {
        const std::string::size_type bar = normalized.find('|');
        if (bar != std::string::npos) {
            normalized.erase(bar);
        }
    }

    if (options.ignore_version) {
        const std::string::size_type dot = normalized.find('.');
        if (dot != std::string::npos) {
            normalized.erase(dot);
        }
    }

    return normalized;
}

} // namespace iso2gene
