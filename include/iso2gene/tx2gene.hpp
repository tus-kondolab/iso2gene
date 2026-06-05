#pragma once

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

#include "iso2gene/id.hpp"
#include "iso2gene/logging.hpp"

namespace iso2gene {

struct Tx2GeneMap {
    std::unordered_map<std::string, std::string> tx_to_gene;
    std::vector<std::string> gene_ids;
    std::unordered_map<std::string, std::size_t> gene_index;
};

Tx2GeneMap read_tx2gene(
    const std::string& path,
    const IdOptions& id_options,
    Logger& logger
);

} // namespace iso2gene
