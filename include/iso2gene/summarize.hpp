#pragma once

#include <string>
#include <vector>

#include "iso2gene/cli.hpp"
#include "iso2gene/logging.hpp"
#include "iso2gene/matrix.hpp"
#include "iso2gene/quant.hpp"
#include "iso2gene/tx2gene.hpp"

namespace iso2gene {

struct GeneMatrices {
    std::vector<std::string> gene_ids;
    std::vector<std::string> sample_names;
    Matrix<double> counts;
    Matrix<double> tpm;
    Matrix<double> length;
    std::size_t total_records = 0;
    std::size_t mapped_records = 0;
    std::size_t unmapped_records = 0;
};

struct TranscriptMatrices {
    std::vector<std::string> transcript_ids;
    std::vector<std::string> gene_ids;
    std::vector<std::string> sample_names;
    Matrix<double> counts;
    Matrix<double> tpm;
    Matrix<double> length;
    std::size_t total_records = 0;
    std::size_t mapped_records = 0;
    std::size_t unmapped_records = 0;
};

GeneMatrices summarize_to_gene(
    const std::vector<SampleInput>& samples,
    const std::vector<std::vector<QuantRecord>>& quantifications,
    const Tx2GeneMap& tx2gene,
    CountMode mode,
    Logger& logger
);

TranscriptMatrices summarize_to_transcripts(
    const std::vector<SampleInput>& samples,
    const std::vector<std::vector<QuantRecord>>& quantifications,
    const Tx2GeneMap& tx2gene,
    TranscriptCountMode mode,
    Logger& logger
);

} // namespace iso2gene
