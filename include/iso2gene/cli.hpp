#pragma once

#include <string>
#include <vector>

#include "iso2gene/id.hpp"

namespace iso2gene {

enum class Command {
    help,
    version,
    counts,
    make_map
};

enum class CountMode {
    simple_sum,
    scaled_tpm,
    length_scaled_tpm
};

struct SampleInput {
    std::string name;
    std::string path;
};

struct Config {
    Command command = Command::help;
    std::string input_type = "kallisto";
    std::string tx2gene_path;
    std::string sample_sheet_path;
    std::string outdir = "out";
    std::string gtf_path;
    std::string map_out_path;
    std::string transcript_id_attr = "transcript_id";
    std::string gene_id_attr = "gene_id";
    CountMode mode = CountMode::length_scaled_tpm;
    IdOptions id_options;
    int precision = 10;
    bool show_help = false;
    std::vector<SampleInput> direct_samples;
};

Config parse_args(int argc, char** argv);

std::string help_text();
std::string version_text();
std::string mode_name(CountMode mode);

} // namespace iso2gene
