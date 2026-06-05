#pragma once

#include <cstddef>
#include <string>

#include "iso2gene/logging.hpp"

namespace iso2gene {

struct GtfMapOptions {
    std::string gtf_path;
    std::string out_path;
    std::string transcript_id_attr = "transcript_id";
    std::string gene_id_attr = "gene_id";
};

struct GtfMapStats {
    std::size_t total_rows = 0;
    std::size_t comment_rows = 0;
    std::size_t empty_rows = 0;
    std::size_t candidate_rows = 0;
    std::size_t missing_transcript_id_rows = 0;
    std::size_t missing_gene_id_rows = 0;
    std::size_t empty_transcript_id_rows = 0;
    std::size_t empty_gene_id_rows = 0;
    std::size_t duplicate_same_gene_rows = 0;
    std::size_t mappings_written = 0;
};

GtfMapStats make_tx2gene_from_gtf(const GtfMapOptions& options, Logger& logger);

} // namespace iso2gene
