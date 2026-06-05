#include "iso2gene/salmon.hpp"

namespace iso2gene {

std::vector<QuantRecord> read_salmon(
    const std::string& path,
    const IdOptions& id_options
) {
    const QuantColumnSpec spec{
        "Salmon quant.sf",
        "Name",
        "Length",
        "EffectiveLength",
        "NumReads",
        "TPM"
    };
    return read_quant_tsv(path, spec, id_options);
}

} // namespace iso2gene
