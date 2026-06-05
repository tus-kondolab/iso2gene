#pragma once

#include <string>
#include <vector>

#include "iso2gene/id.hpp"
#include "iso2gene/quant.hpp"

namespace iso2gene {

std::vector<QuantRecord> read_quantification_file(
    const std::string& input_type,
    const std::string& path,
    const IdOptions& id_options
);

} // namespace iso2gene
