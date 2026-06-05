#include "iso2gene/gtf.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "iso2gene/error.hpp"
#include "iso2gene/text_reader.hpp"
#include "iso2gene/tsv.hpp"

namespace iso2gene {

namespace {

struct AttributeValues {
    bool has_transcript_id = false;
    bool has_gene_id = false;
    bool saw_empty_transcript_id = false;
    bool saw_empty_gene_id = false;
    std::string transcript_id;
    std::string gene_id;
};

struct MappingEntry {
    std::string transcript_id;
    std::string gene_id;
};

std::string trim_ascii(const std::string& value) {
    std::string::size_type begin = 0;
    while (begin < value.size() && std::isspace(static_cast<unsigned char>(value[begin]))) {
        ++begin;
    }

    std::string::size_type end = value.size();
    while (end > begin && std::isspace(static_cast<unsigned char>(value[end - 1]))) {
        --end;
    }

    return value.substr(begin, end - begin);
}

bool starts_with_token(const std::string& line, const std::string& token) {
    if (line.size() < token.size() || line.compare(0, token.size(), token) != 0) {
        return false;
    }
    return line.size() == token.size()
        || std::isspace(static_cast<unsigned char>(line[token.size()]));
}

std::string unquote_value(const std::string& raw_value) {
    std::string value = trim_ascii(raw_value);
    if (value.size() >= 2 && value.front() == '"' && value.back() == '"') {
        value = value.substr(1, value.size() - 2);
    }
    return value;
}

void assign_required_attribute(
    const std::string& path,
    std::size_t line_number,
    const std::string& attr_name,
    const std::string& value,
    bool& has_value,
    bool& saw_empty,
    std::string& current_value
) {
    if (value.empty()) {
        saw_empty = true;
        return;
    }
    if (has_value && current_value != value) {
        std::ostringstream msg;
        msg << path << ":" << line_number << ": attribute '" << attr_name
            << "' appears multiple times with different values ('" << current_value
            << "' and '" << value << "')";
        throw Iso2GeneError(ExitCode::parse_error, msg.str());
    }
    has_value = true;
    current_value = value;
}

AttributeValues parse_gtf_attributes(
    const std::string& path,
    std::size_t line_number,
    const std::string& attributes,
    const GtfMapOptions& options
) {
    AttributeValues values;
    std::string::size_type start = 0;
    while (start <= attributes.size()) {
        const std::string::size_type semicolon = attributes.find(';', start);
        const std::string raw_item = semicolon == std::string::npos
            ? attributes.substr(start)
            : attributes.substr(start, semicolon - start);
        const std::string item = trim_ascii(raw_item);

        if (!item.empty()) {
            std::string::size_type split = 0;
            while (split < item.size()
                   && !std::isspace(static_cast<unsigned char>(item[split]))) {
                ++split;
            }
            const std::string key = item.substr(0, split);
            const std::string value = split < item.size()
                ? unquote_value(item.substr(split))
                : std::string();

            if (key == options.transcript_id_attr) {
                assign_required_attribute(
                    path,
                    line_number,
                    key,
                    value,
                    values.has_transcript_id,
                    values.saw_empty_transcript_id,
                    values.transcript_id
                );
            } else if (key == options.gene_id_attr) {
                assign_required_attribute(
                    path,
                    line_number,
                    key,
                    value,
                    values.has_gene_id,
                    values.saw_empty_gene_id,
                    values.gene_id
                );
            }
        }

        if (semicolon == std::string::npos) {
            break;
        }
        start = semicolon + 1;
    }
    return values;
}

void write_tx2gene_tsv(
    const std::string& path,
    const std::vector<MappingEntry>& mappings
) {
    std::ofstream out(path);
    if (!out) {
        throw Iso2GeneError(ExitCode::io_error, "failed to open tx2gene output for writing: " + path);
    }

    out << "transcript_id\tgene_id\n";
    for (const MappingEntry& entry : mappings) {
        out << entry.transcript_id << '\t' << entry.gene_id << '\n';
    }
}

void log_gtf_stats(const GtfMapStats& stats, Logger& logger) {
    logger.info("parsed " + std::to_string(stats.total_rows) + " GTF rows");
    logger.info("wrote " + std::to_string(stats.mappings_written) + " transcript-to-gene mappings");
    if (stats.missing_transcript_id_rows > 0) {
        logger.warn(
            "ignored " + std::to_string(stats.missing_transcript_id_rows)
            + " rows without transcript_id"
        );
    }
    if (stats.missing_gene_id_rows > 0) {
        logger.warn(
            "ignored " + std::to_string(stats.missing_gene_id_rows)
            + " rows without gene_id"
        );
    }
    if (stats.empty_transcript_id_rows > 0) {
        logger.warn(
            "ignored " + std::to_string(stats.empty_transcript_id_rows)
            + " rows with empty transcript_id"
        );
    }
    if (stats.empty_gene_id_rows > 0) {
        logger.warn(
            "ignored " + std::to_string(stats.empty_gene_id_rows)
            + " rows with empty gene_id"
        );
    }
    if (stats.duplicate_same_gene_rows > 0) {
        logger.warn(
            "ignored " + std::to_string(stats.duplicate_same_gene_rows)
            + " duplicate transcript-to-gene rows"
        );
    }
}

} // namespace

GtfMapStats make_tx2gene_from_gtf(const GtfMapOptions& options, Logger& logger) {
    if (options.gtf_path.empty()) {
        throw Iso2GeneError(ExitCode::input_error, "missing required option: --gtf");
    }
    if (options.out_path.empty()) {
        throw Iso2GeneError(ExitCode::input_error, "missing required option: --out");
    }
    if (options.transcript_id_attr.empty()) {
        throw Iso2GeneError(ExitCode::input_error, "--transcript-id-attr must be non-empty");
    }
    if (options.gene_id_attr.empty()) {
        throw Iso2GeneError(ExitCode::input_error, "--gene-id-attr must be non-empty");
    }

    TextReader input(options.gtf_path);

    GtfMapStats stats;
    std::vector<MappingEntry> mappings;
    std::unordered_map<std::string, std::string> tx_to_gene;

    std::string line;
    while (input.read_line(line)) {
        ++stats.total_rows;

        const std::string trimmed_line = trim_ascii(line);
        if (trimmed_line.empty()) {
            ++stats.empty_rows;
            continue;
        }
        if (trimmed_line[0] == '#') {
            ++stats.comment_rows;
            continue;
        }
        if (starts_with_token(trimmed_line, "track") || starts_with_token(trimmed_line, "browser")) {
            std::ostringstream msg;
            msg << options.gtf_path << ":" << stats.total_rows
                << ": unsupported UCSC track/browser directive; expected a 9-column GTF feature line";
            throw Iso2GeneError(ExitCode::parse_error, msg.str());
        }

        const std::vector<std::string> fields = split_tsv_line(line);
        if (fields.size() != 9) {
            std::ostringstream msg;
            msg << options.gtf_path << ":" << stats.total_rows
                << ": expected 9 GTF columns, found " << fields.size()
                << "; check that the input is GTF, not GFF3";
            throw Iso2GeneError(ExitCode::parse_error, msg.str());
        }

        const AttributeValues attrs =
            parse_gtf_attributes(options.gtf_path, stats.total_rows, fields[8], options);

        if (!attrs.has_transcript_id) {
            if (attrs.saw_empty_transcript_id) {
                ++stats.empty_transcript_id_rows;
            } else {
                ++stats.missing_transcript_id_rows;
            }
            continue;
        }
        if (!attrs.has_gene_id) {
            if (attrs.saw_empty_gene_id) {
                ++stats.empty_gene_id_rows;
            } else {
                ++stats.missing_gene_id_rows;
            }
            continue;
        }

        ++stats.candidate_rows;
        const auto existing = tx_to_gene.find(attrs.transcript_id);
        if (existing != tx_to_gene.end()) {
            if (existing->second != attrs.gene_id) {
                std::ostringstream msg;
                msg << options.gtf_path << ":" << stats.total_rows
                    << ": transcript '" << attrs.transcript_id
                    << "' maps to multiple genes ('" << existing->second
                    << "' and '" << attrs.gene_id << "')";
                throw Iso2GeneError(ExitCode::parse_error, msg.str());
            }
            ++stats.duplicate_same_gene_rows;
            continue;
        }

        tx_to_gene.emplace(attrs.transcript_id, attrs.gene_id);
        mappings.push_back(MappingEntry{attrs.transcript_id, attrs.gene_id});
    }

    if (mappings.empty()) {
        throw Iso2GeneError(
            ExitCode::parse_error,
            options.gtf_path
            + ": no transcript-to-gene mappings found; check that the input is GTF "
              "(GFF3 is not supported) and contains transcript_id/gene_id attributes"
        );
    }

    write_tx2gene_tsv(options.out_path, mappings);
    stats.mappings_written = mappings.size();
    log_gtf_stats(stats, logger);
    return stats;
}

} // namespace iso2gene
