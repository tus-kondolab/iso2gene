#include "iso2gene/sample_sheet.hpp"

#include <sstream>
#include <unordered_set>

#include "iso2gene/error.hpp"
#include "iso2gene/tsv.hpp"

namespace iso2gene {

std::vector<SampleInput> read_sample_sheet(const std::string& path) {
    TsvReader reader(path);
    std::vector<std::string> header;
    if (!reader.read_row(header)) {
        throw Iso2GeneError(ExitCode::parse_error, path + ": empty sample sheet");
    }

    const int sample_col = find_column(header, "sample", path);
    const int path_col = find_column(header, "path", path);
    const std::size_t min_cols = header.size();

    std::vector<SampleInput> samples;
    std::unordered_set<std::string> seen_names;
    std::vector<std::string> fields;
    while (reader.read_row(fields)) {
        if (fields.size() != min_cols) {
            std::ostringstream msg;
            msg << path << ":" << reader.line_number()
                << ": expected " << min_cols << " columns, found " << fields.size();
            throw Iso2GeneError(ExitCode::parse_error, msg.str());
        }

        const std::string& sample = fields[static_cast<std::size_t>(sample_col)];
        const std::string& quant_path = fields[static_cast<std::size_t>(path_col)];
        if (sample.empty()) {
            std::ostringstream msg;
            msg << path << ":" << reader.line_number() << ": empty sample name";
            throw Iso2GeneError(ExitCode::parse_error, msg.str());
        }
        if (quant_path.empty()) {
            std::ostringstream msg;
            msg << path << ":" << reader.line_number() << ": empty path";
            throw Iso2GeneError(ExitCode::parse_error, msg.str());
        }
        if (!seen_names.insert(sample).second) {
            std::ostringstream msg;
            msg << path << ":" << reader.line_number()
                << ": duplicate sample name: " << sample;
            throw Iso2GeneError(ExitCode::parse_error, msg.str());
        }
        samples.push_back(SampleInput{sample, quant_path});
    }

    if (samples.empty()) {
        throw Iso2GeneError(ExitCode::parse_error, path + ": sample sheet has no samples");
    }
    return samples;
}

} // namespace iso2gene
