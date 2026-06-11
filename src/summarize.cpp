#include "iso2gene/summarize.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "iso2gene/error.hpp"

namespace iso2gene {

namespace {

struct LengthAccumulator {
    double weighted_sum = 0.0;
    double tpm_sum = 0.0;
};

struct TranscriptLengthAccumulator {
    std::string gene_id;
    double sum = 0.0;
    std::size_t count = 0;
};

struct MeanAccumulator {
    double sum = 0.0;
    std::size_t count = 0;
};

bool is_missing_length(double value) {
    return std::isnan(value);
}

double missing_length() {
    return std::numeric_limits<double>::quiet_NaN();
}

double row_mean(const Matrix<double>& matrix, std::size_t row) {
    double sum = 0.0;
    for (std::size_t col = 0; col < matrix.cols(); ++col) {
        sum += matrix(row, col);
    }
    return matrix.cols() > 0 ? sum / static_cast<double>(matrix.cols()) : 0.0;
}

double geometric_mean_non_missing(const Matrix<double>& matrix, std::size_t row) {
    double log_sum = 0.0;
    std::size_t count = 0;
    bool has_zero = false;

    for (std::size_t col = 0; col < matrix.cols(); ++col) {
        const double value = matrix(row, col);
        if (is_missing_length(value)) {
            continue;
        }
        if (value <= 0.0) {
            has_zero = true;
            ++count;
            continue;
        }
        log_sum += std::log(value);
        ++count;
    }

    if (count == 0 || has_zero) {
        return 0.0;
    }
    return std::exp(log_sum / static_cast<double>(count));
}

void replace_missing_lengths(Matrix<double>& length, const std::vector<double>& average_gene_length) {
    for (std::size_t row = 0; row < length.rows(); ++row) {
        bool any_missing = false;
        bool all_missing = true;
        for (std::size_t col = 0; col < length.cols(); ++col) {
            if (is_missing_length(length(row, col))) {
                any_missing = true;
            } else {
                all_missing = false;
            }
        }

        if (!any_missing) {
            continue;
        }

        const double replacement = all_missing
            ? average_gene_length[row]
            : geometric_mean_non_missing(length, row);

        for (std::size_t col = 0; col < length.cols(); ++col) {
            if (is_missing_length(length(row, col))) {
                length(row, col) = replacement;
            }
        }
    }
}

double column_sum(const Matrix<double>& matrix, std::size_t col) {
    double total = 0.0;
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        total += matrix(row, col);
    }
    return total;
}

double median(std::vector<double> values) {
    if (values.empty()) {
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    const std::size_t mid = values.size() / 2;
    if (values.size() % 2 == 1) {
        return values[mid];
    }
    return (values[mid - 1] + values[mid]) / 2.0;
}

void counts_from_abundance(GeneMatrices& output, CountMode mode) {
    if (mode == CountMode::simple_sum) {
        return;
    }

    std::vector<double> mean_gene_lengths(output.gene_ids.size(), 0.0);
    if (mode == CountMode::length_scaled_tpm) {
        for (std::size_t gene = 0; gene < output.gene_ids.size(); ++gene) {
            mean_gene_lengths[gene] = row_mean(output.length, gene);
        }
    }

    for (std::size_t sample = 0; sample < output.sample_names.size(); ++sample) {
        const double counts_sum = column_sum(output.counts, sample);

        double new_sum = 0.0;
        for (std::size_t gene = 0; gene < output.gene_ids.size(); ++gene) {
            const double new_count = mode == CountMode::scaled_tpm
                ? output.tpm(gene, sample)
                : output.tpm(gene, sample) * mean_gene_lengths[gene];
            new_sum += new_count;
        }

        const double scale = new_sum > 0.0 ? counts_sum / new_sum : 0.0;
        for (std::size_t gene = 0; gene < output.gene_ids.size(); ++gene) {
            const double new_count = mode == CountMode::scaled_tpm
                ? output.tpm(gene, sample)
                : output.tpm(gene, sample) * mean_gene_lengths[gene];
            output.counts(gene, sample) = new_count * scale;
        }
    }
}

} // namespace

GeneMatrices summarize_to_gene(
    const std::vector<SampleInput>& samples,
    const std::vector<std::vector<QuantRecord>>& quantifications,
    const Tx2GeneMap& tx2gene,
    CountMode mode,
    Logger& logger
) {
    if (samples.empty()) {
        throw Iso2GeneError(ExitCode::input_error, "no samples provided");
    }
    if (samples.size() != quantifications.size()) {
        throw Iso2GeneError(
            ExitCode::internal_error,
            "sample list and quantification list sizes differ"
        );
    }

    std::set<std::string> observed_genes;
    std::unordered_map<std::string, TranscriptLengthAccumulator> transcript_lengths;
    for (const std::vector<QuantRecord>& records : quantifications) {
        for (const QuantRecord& record : records) {
            const auto gene_it = tx2gene.tx_to_gene.find(record.transcript_id);
            if (gene_it == tx2gene.tx_to_gene.end()) {
                continue;
            }
            observed_genes.insert(gene_it->second);
            TranscriptLengthAccumulator& acc = transcript_lengths[record.transcript_id];
            acc.gene_id = gene_it->second;
            acc.sum += record.eff_length;
            ++acc.count;
        }
    }

    if (observed_genes.empty()) {
        throw Iso2GeneError(
            ExitCode::input_error,
            "no quantification transcript records matched the tx2gene map"
        );
    }

    GeneMatrices output;
    output.gene_ids.assign(observed_genes.begin(), observed_genes.end());
    std::unordered_map<std::string, std::size_t> observed_gene_index;
    for (std::size_t gene = 0; gene < output.gene_ids.size(); ++gene) {
        observed_gene_index.emplace(output.gene_ids[gene], gene);
    }

    output.sample_names.reserve(samples.size());
    for (const SampleInput& sample : samples) {
        output.sample_names.push_back(sample.name);
    }
    output.counts = Matrix<double>(output.gene_ids.size(), samples.size(), 0.0);
    output.tpm = Matrix<double>(output.gene_ids.size(), samples.size(), 0.0);
    output.length = Matrix<double>(output.gene_ids.size(), samples.size(), missing_length());

    std::vector<MeanAccumulator> average_gene_length_acc(output.gene_ids.size());
    for (const auto& item : transcript_lengths) {
        const TranscriptLengthAccumulator& tx_len = item.second;
        if (tx_len.count == 0) {
            continue;
        }
        const auto gene_index_it = observed_gene_index.find(tx_len.gene_id);
        if (gene_index_it == observed_gene_index.end()) {
            continue;
        }
        MeanAccumulator& gene_acc = average_gene_length_acc[gene_index_it->second];
        gene_acc.sum += tx_len.sum / static_cast<double>(tx_len.count);
        ++gene_acc.count;
    }

    std::vector<double> average_gene_length(output.gene_ids.size(), 0.0);
    for (std::size_t gene = 0; gene < output.gene_ids.size(); ++gene) {
        const MeanAccumulator& acc = average_gene_length_acc[gene];
        average_gene_length[gene] =
            acc.count > 0 ? acc.sum / static_cast<double>(acc.count) : 0.0;
    }

    for (std::size_t sample_index = 0; sample_index < samples.size(); ++sample_index) {
        const std::vector<QuantRecord>& records = quantifications[sample_index];

        std::vector<LengthAccumulator> length_acc(output.gene_ids.size());

        for (const QuantRecord& record : records) {
            ++output.total_records;

            const auto gene_it = tx2gene.tx_to_gene.find(record.transcript_id);
            if (gene_it == tx2gene.tx_to_gene.end()) {
                ++output.unmapped_records;
                continue;
            }

            const auto gene_index_it = observed_gene_index.find(gene_it->second);
            if (gene_index_it == observed_gene_index.end()) {
                throw Iso2GeneError(
                    ExitCode::internal_error,
                    "gene index missing for gene: " + gene_it->second
                );
            }
            const std::size_t gene_index = gene_index_it->second;
            ++output.mapped_records;

            output.tpm(gene_index, sample_index) += record.tpm;
            output.counts(gene_index, sample_index) += record.est_counts;

            LengthAccumulator& len = length_acc[gene_index];
            len.weighted_sum += record.tpm * record.eff_length;
            len.tpm_sum += record.tpm;
        }

        for (std::size_t gene_index = 0; gene_index < output.gene_ids.size(); ++gene_index) {
            const LengthAccumulator& len = length_acc[gene_index];
            if (len.tpm_sum > 0.0) {
                output.length(gene_index, sample_index) = len.weighted_sum / len.tpm_sum;
            } else {
                output.length(gene_index, sample_index) = missing_length();
            }
        }
    }

    if (output.mapped_records == 0) {
        throw Iso2GeneError(
            ExitCode::input_error,
            "no quantification transcript records matched the tx2gene map"
        );
    }
    if (output.unmapped_records > 0) {
        logger.warn(
            std::to_string(output.unmapped_records)
            + " transcript records were not found in tx2gene and were ignored"
        );
    }

    replace_missing_lengths(output.length, average_gene_length);
    counts_from_abundance(output, mode);

    return output;
}

TranscriptMatrices summarize_to_transcripts(
    const std::vector<SampleInput>& samples,
    const std::vector<std::vector<QuantRecord>>& quantifications,
    const Tx2GeneMap& tx2gene,
    TranscriptCountMode mode,
    Logger& logger
) {
    if (samples.empty()) {
        throw Iso2GeneError(ExitCode::input_error, "no samples provided");
    }
    if (samples.size() != quantifications.size()) {
        throw Iso2GeneError(
            ExitCode::internal_error,
            "sample list and quantification list sizes differ"
        );
    }
    if (mode != TranscriptCountMode::dtu_scaled_tpm) {
        throw Iso2GeneError(ExitCode::internal_error, "unsupported transcript count mode");
    }

    TranscriptMatrices output;
    output.sample_names.reserve(samples.size());
    for (const SampleInput& sample : samples) {
        output.sample_names.push_back(sample.name);
    }

    std::unordered_map<std::string, std::size_t> transcript_index;
    for (const std::vector<QuantRecord>& records : quantifications) {
        for (const QuantRecord& record : records) {
            ++output.total_records;

            const auto gene_it = tx2gene.tx_to_gene.find(record.transcript_id);
            if (gene_it == tx2gene.tx_to_gene.end()) {
                ++output.unmapped_records;
                continue;
            }
            ++output.mapped_records;

            if (transcript_index.find(record.transcript_id) == transcript_index.end()) {
                const std::size_t index = output.transcript_ids.size();
                transcript_index.emplace(record.transcript_id, index);
                output.transcript_ids.push_back(record.transcript_id);
                output.gene_ids.push_back(gene_it->second);
            }
        }
    }

    if (output.transcript_ids.empty() || output.mapped_records == 0) {
        throw Iso2GeneError(
            ExitCode::input_error,
            "no quantification transcript records matched the tx2gene map"
        );
    }

    const std::size_t rows = output.transcript_ids.size();
    const std::size_t cols = samples.size();
    output.counts = Matrix<double>(rows, cols, 0.0);
    output.tpm = Matrix<double>(rows, cols, 0.0);
    output.length = Matrix<double>(rows, cols, 0.0);
    std::vector<unsigned char> observed(rows * cols, 0);
    std::vector<double> library_sizes(cols, 0.0);

    for (std::size_t sample_index = 0; sample_index < samples.size(); ++sample_index) {
        const std::vector<QuantRecord>& records = quantifications[sample_index];
        for (const QuantRecord& record : records) {
            const auto index_it = transcript_index.find(record.transcript_id);
            if (index_it == transcript_index.end()) {
                continue;
            }

            const std::size_t row = index_it->second;
            output.counts(row, sample_index) = record.est_counts;
            output.tpm(row, sample_index) = record.tpm;
            output.length(row, sample_index) = record.eff_length;
            observed[row * cols + sample_index] = 1;
            library_sizes[sample_index] += record.est_counts;
        }
    }

    std::unordered_map<std::string, std::vector<double>> gene_mean_lengths;
    for (std::size_t row = 0; row < rows; ++row) {
        double sum = 0.0;
        std::size_t count = 0;
        for (std::size_t col = 0; col < cols; ++col) {
            if (observed[row * cols + col]) {
                sum += output.length(row, col);
                ++count;
            }
        }
        const double mean_length = count > 0 ? sum / static_cast<double>(count) : 0.0;
        gene_mean_lengths[output.gene_ids[row]].push_back(mean_length);
    }

    std::unordered_map<std::string, double> gene_median_lengths;
    for (const auto& item : gene_mean_lengths) {
        gene_median_lengths.emplace(item.first, median(item.second));
    }

    for (std::size_t sample_index = 0; sample_index < cols; ++sample_index) {
        double raw_sum = 0.0;
        for (std::size_t row = 0; row < rows; ++row) {
            const double median_length = gene_median_lengths[output.gene_ids[row]];
            raw_sum += output.tpm(row, sample_index) * median_length;
        }

        const double scale = raw_sum > 0.0 ? library_sizes[sample_index] / raw_sum : 0.0;
        for (std::size_t row = 0; row < rows; ++row) {
            const double median_length = gene_median_lengths[output.gene_ids[row]];
            output.counts(row, sample_index) =
                output.tpm(row, sample_index) * median_length * scale;
        }
    }

    if (output.unmapped_records > 0) {
        logger.warn(
            std::to_string(output.unmapped_records)
            + " transcript records were not found in tx2gene and were ignored"
        );
    }

    return output;
}

} // namespace iso2gene
