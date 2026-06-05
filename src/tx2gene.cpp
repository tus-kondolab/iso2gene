#include "iso2gene/tx2gene.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>

#include "iso2gene/error.hpp"
#include "iso2gene/tsv.hpp"

namespace iso2gene {

namespace {

std::string lower_ascii(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

bool looks_like_header(const std::vector<std::string>& row) {
    if (row.size() < 2) {
        return false;
    }
    const std::string first = lower_ascii(row[0]);
    const std::string second = lower_ascii(row[1]);
    const bool first_ok =
        first == "transcript_id" || first == "transcript" || first == "txname" || first == "tx_id";
    const bool second_ok =
        second == "gene_id" || second == "gene" || second == "geneid";
    return first_ok && second_ok;
}

void add_gene_if_new(Tx2GeneMap& map, const std::string& gene_id) {
    if (map.gene_index.find(gene_id) != map.gene_index.end()) {
        return;
    }
    const std::size_t index = map.gene_ids.size();
    map.gene_ids.push_back(gene_id);
    map.gene_index.emplace(gene_id, index);
}

void read_mapping_row(
    Tx2GeneMap& map,
    Logger& logger,
    const std::string& path,
    std::size_t line_number,
    const std::vector<std::string>& fields,
    const IdOptions& id_options,
    std::size_t& duplicate_same_gene
) {
    if (fields.size() < 2) {
        std::ostringstream msg;
        msg << path << ":" << line_number << ": expected at least 2 columns";
        throw Iso2GeneError(ExitCode::parse_error, msg.str());
    }

    const std::string transcript_id = normalize_transcript_id(fields[0], id_options);
    const std::string& gene_id = fields[1];
    if (transcript_id.empty()) {
        std::ostringstream msg;
        msg << path << ":" << line_number << ": empty transcript_id";
        throw Iso2GeneError(ExitCode::parse_error, msg.str());
    }
    if (gene_id.empty()) {
        std::ostringstream msg;
        msg << path << ":" << line_number << ": empty gene_id";
        throw Iso2GeneError(ExitCode::parse_error, msg.str());
    }

    const auto existing = map.tx_to_gene.find(transcript_id);
    if (existing != map.tx_to_gene.end()) {
        if (existing->second != gene_id) {
            std::ostringstream msg;
            msg << path << ":" << line_number << ": transcript '" << transcript_id
                << "' maps to multiple genes ('" << existing->second << "' and '"
                << gene_id << "')";
            throw Iso2GeneError(ExitCode::parse_error, msg.str());
        }
        ++duplicate_same_gene;
        return;
    }

    map.tx_to_gene.emplace(transcript_id, gene_id);
    add_gene_if_new(map, gene_id);

    (void)logger;
}

} // namespace

Tx2GeneMap read_tx2gene(
    const std::string& path,
    const IdOptions& id_options,
    Logger& logger
) {
    TsvReader reader(path);
    Tx2GeneMap map;
    std::size_t duplicate_same_gene = 0;

    std::vector<std::string> fields;
    if (!reader.read_row(fields)) {
        throw Iso2GeneError(ExitCode::parse_error, path + ": empty tx2gene file");
    }

    if (!looks_like_header(fields)) {
        read_mapping_row(
            map,
            logger,
            path,
            reader.line_number(),
            fields,
            id_options,
            duplicate_same_gene
        );
    }

    while (reader.read_row(fields)) {
        read_mapping_row(
            map,
            logger,
            path,
            reader.line_number(),
            fields,
            id_options,
            duplicate_same_gene
        );
    }

    if (map.tx_to_gene.empty()) {
        throw Iso2GeneError(ExitCode::parse_error, path + ": no transcript-to-gene mappings found");
    }
    if (duplicate_same_gene > 0) {
        logger.warn(
            std::to_string(duplicate_same_gene)
            + " duplicate tx2gene rows mapped the same transcript to the same gene and were ignored"
        );
    }

    return map;
}

} // namespace iso2gene
