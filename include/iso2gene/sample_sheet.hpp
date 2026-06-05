#pragma once

#include <string>
#include <vector>

#include "iso2gene/cli.hpp"

namespace iso2gene {

std::vector<SampleInput> read_sample_sheet(const std::string& path);

} // namespace iso2gene
