#pragma once

#include <string>

namespace iso2gene {

struct IdOptions {
    bool ignore_version = false;
    bool ignore_after_bar = false;
};

std::string normalize_transcript_id(const std::string& id, const IdOptions& options);

} // namespace iso2gene
