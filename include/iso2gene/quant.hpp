#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "iso2gene/id.hpp"

namespace iso2gene {

struct QuantRecord {
    std::string transcript_id;
    double length = 0.0;
    double eff_length = 0.0;
    double est_counts = 0.0;
    double tpm = 0.0;
};

struct QuantColumnSpec {
    const char* format_name;
    const char* transcript_id;
    const char* length;
    const char* eff_length;
    const char* est_counts;
    const char* tpm;
};

std::vector<QuantRecord> read_quant_tsv(
    const std::string& path,
    const QuantColumnSpec& spec,
    const IdOptions& id_options
);

void require_non_negative_quant_value(
    double value,
    const std::string& path,
    std::size_t line_number,
    const std::string& column_name
);

} // namespace iso2gene
