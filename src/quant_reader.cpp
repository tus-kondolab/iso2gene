#include "iso2gene/quant_reader.hpp"

#include "iso2gene/error.hpp"
#include "iso2gene/kallisto.hpp"
#include "iso2gene/rsem.hpp"
#include "iso2gene/salmon.hpp"

namespace iso2gene {

std::vector<QuantRecord> read_quantification_file(
    const std::string& input_type,
    const std::string& path,
    const IdOptions& id_options
) {
    if (input_type == "kallisto") {
        return read_kallisto(path, id_options);
    }
    if (input_type == "salmon") {
        return read_salmon(path, id_options);
    }
    if (input_type == "rsem") {
        return read_rsem(path, id_options);
    }
    throw Iso2GeneError(
        ExitCode::input_error,
        "unsupported --type '" + input_type + "'"
    );
}

} // namespace iso2gene
