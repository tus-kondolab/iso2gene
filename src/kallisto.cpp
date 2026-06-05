#include "iso2gene/kallisto.hpp"

namespace iso2gene {

std::vector<QuantRecord> read_kallisto(
    const std::string& path,
    const IdOptions& id_options
) {
    const QuantColumnSpec spec{
        "kallisto abundance",
        "target_id",
        "length",
        "eff_length",
        "est_counts",
        "tpm"
    };
    return read_quant_tsv(path, spec, id_options);
}

} // namespace iso2gene
