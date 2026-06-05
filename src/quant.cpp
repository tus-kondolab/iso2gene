#include "iso2gene/quant.hpp"

#include <cstddef>
#include <sstream>
#include <unordered_set>

#include "iso2gene/error.hpp"
#include "iso2gene/tsv.hpp"

namespace iso2gene {

std::vector<QuantRecord> read_quant_tsv(
    const std::string& path,
    const QuantColumnSpec& spec,
    const IdOptions& id_options
) {
    TsvReader reader(path);
    std::vector<std::string> header;
    if (!reader.read_row(header)) {
        throw Iso2GeneError(
            ExitCode::parse_error,
            path + ": empty " + spec.format_name + " file"
        );
    }

    const int transcript_col = find_column(header, spec.transcript_id, path);
    const int length_col = find_column(header, spec.length, path);
    const int eff_length_col = find_column(header, spec.eff_length, path);
    const int est_counts_col = find_column(header, spec.est_counts, path);
    const int tpm_col = find_column(header, spec.tpm, path);
    const std::size_t expected_cols = header.size();

    std::vector<QuantRecord> records;
    std::unordered_set<std::string> seen_transcripts;
    std::vector<std::string> fields;
    while (reader.read_row(fields)) {
        if (fields.size() != expected_cols) {
            std::ostringstream msg;
            msg << path << ":" << reader.line_number()
                << ": expected " << expected_cols << " columns, found " << fields.size();
            throw Iso2GeneError(ExitCode::parse_error, msg.str());
        }

        QuantRecord record;
        record.transcript_id = normalize_transcript_id(
            fields[static_cast<std::size_t>(transcript_col)],
            id_options
        );
        if (record.transcript_id.empty()) {
            std::ostringstream msg;
            msg << path << ":" << reader.line_number()
                << ": empty " << spec.transcript_id;
            throw Iso2GeneError(ExitCode::parse_error, msg.str());
        }
        if (!seen_transcripts.insert(record.transcript_id).second) {
            std::ostringstream msg;
            msg << path << ":" << reader.line_number()
                << ": duplicate " << spec.transcript_id
                << " after normalization: " << record.transcript_id;
            throw Iso2GeneError(ExitCode::parse_error, msg.str());
        }

        record.length = parse_double_strict(
            fields[static_cast<std::size_t>(length_col)],
            path,
            reader.line_number(),
            spec.length
        );
        record.eff_length = parse_double_strict(
            fields[static_cast<std::size_t>(eff_length_col)],
            path,
            reader.line_number(),
            spec.eff_length
        );
        record.est_counts = parse_double_strict(
            fields[static_cast<std::size_t>(est_counts_col)],
            path,
            reader.line_number(),
            spec.est_counts
        );
        record.tpm = parse_double_strict(
            fields[static_cast<std::size_t>(tpm_col)],
            path,
            reader.line_number(),
            spec.tpm
        );

        require_non_negative_quant_value(record.length, path, reader.line_number(), spec.length);
        require_non_negative_quant_value(record.eff_length, path, reader.line_number(), spec.eff_length);
        require_non_negative_quant_value(record.est_counts, path, reader.line_number(), spec.est_counts);
        require_non_negative_quant_value(record.tpm, path, reader.line_number(), spec.tpm);

        records.push_back(record);
    }

    if (records.empty()) {
        throw Iso2GeneError(
            ExitCode::parse_error,
            path + ": no " + spec.format_name + " records found"
        );
    }

    return records;
}

void require_non_negative_quant_value(
    double value,
    const std::string& path,
    std::size_t line_number,
    const std::string& column_name
) {
    if (value < 0.0) {
        std::ostringstream msg;
        msg << path << ":" << line_number << ": negative value in '" << column_name << "'";
        throw Iso2GeneError(ExitCode::parse_error, msg.str());
    }
}

} // namespace iso2gene
