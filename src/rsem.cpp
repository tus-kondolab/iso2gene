#include "iso2gene/rsem.hpp"

namespace iso2gene {

std::vector<QuantRecord> read_rsem(
    const std::string& path,
    const IdOptions& id_options
) {
    const QuantColumnSpec spec{
        "RSEM isoforms.results",
        "transcript_id",
        "length",
        "effective_length",
        "expected_count",
        "TPM"
    };
    std::vector<QuantRecord> records = read_quant_tsv(path, spec, id_options);
    for (QuantRecord& record : records) {
        // tximport 1.38.2 applies this RSEM-specific clamp before
        // gene-level summarization.
        if (record.eff_length < 1.0) {
            record.eff_length = 1.0;
        }
    }
    return records;
}

} // namespace iso2gene
