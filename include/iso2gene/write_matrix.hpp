#pragma once

#include <string>
#include <vector>

#include "iso2gene/cli.hpp"
#include "iso2gene/logging.hpp"
#include "iso2gene/matrix.hpp"
#include "iso2gene/summarize.hpp"

namespace iso2gene {

void write_matrix_tsv(
    const std::string& path,
    const std::string& row_name,
    const std::vector<std::string>& row_ids,
    const std::vector<std::string>& col_names,
    const Matrix<double>& matrix,
    int precision
);

void write_outputs(
    const Config& config,
    const GeneMatrices& matrices,
    const Logger& logger
);

} // namespace iso2gene
