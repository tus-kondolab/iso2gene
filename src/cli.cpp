#include "iso2gene/cli.hpp"

#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

#include "iso2gene/error.hpp"
#include "iso2gene/version.hpp"

namespace iso2gene {

namespace {

bool is_option(const std::string& arg) {
    return arg.size() >= 2 && arg[0] == '-' && arg[1] == '-';
}

std::string require_value(int argc, char** argv, int& index, const std::string& option) {
    if (index + 1 >= argc) {
        throw Iso2GeneError(ExitCode::input_error, "missing value for option: " + option);
    }
    ++index;
    return argv[index];
}

CountMode parse_mode(const std::string& value) {
    if (value == "simple-sum") {
        return CountMode::simple_sum;
    }
    if (value == "scaled-tpm") {
        return CountMode::scaled_tpm;
    }
    if (value == "length-scaled-tpm") {
        return CountMode::length_scaled_tpm;
    }
    if (value == "dtu-scaled-tpm" || value == "dtu_scaled_tpm") {
        throw Iso2GeneError(
            ExitCode::input_error,
            "mode '" + value + "' is transcript-level only; use `iso2gene txout --mode dtu-scaled-tpm`"
        );
    }
    throw Iso2GeneError(
        ExitCode::input_error,
        "unsupported mode '" + value + "'; expected simple-sum, scaled-tpm, or length-scaled-tpm"
    );
}

TranscriptCountMode parse_transcript_mode(const std::string& value) {
    if (value == "dtu-scaled-tpm" || value == "dtu_scaled_tpm") {
        return TranscriptCountMode::dtu_scaled_tpm;
    }
    throw Iso2GeneError(
        ExitCode::input_error,
        "unsupported txout mode '" + value + "'; expected dtu-scaled-tpm"
    );
}

SampleInput parse_direct_sample(const std::string& arg) {
    const std::string::size_type eq = arg.find('=');
    if (eq == std::string::npos || eq == 0 || eq == arg.size() - 1) {
        throw Iso2GeneError(
            ExitCode::input_error,
            "direct input must have the form sample=path: " + arg
        );
    }
    return SampleInput{arg.substr(0, eq), arg.substr(eq + 1)};
}

int parse_precision(const std::string& value) {
    char* end = nullptr;
    const long parsed = std::strtol(value.c_str(), &end, 10);
    if (end == value.c_str() || *end != '\0' || parsed < 1 || parsed > 17) {
        throw Iso2GeneError(
            ExitCode::input_error,
            "invalid --precision value; expected integer from 1 to 17"
        );
    }
    return static_cast<int>(parsed);
}

void validate_quantification_config(const Config& config) {
    if (
        config.input_type != "kallisto"
        && config.input_type != "salmon"
        && config.input_type != "rsem"
    ) {
        throw Iso2GeneError(
            ExitCode::input_error,
            "unsupported --type '" + config.input_type + "'; expected kallisto, salmon, or rsem"
        );
    }
    if (config.tx2gene_path.empty()) {
        throw Iso2GeneError(ExitCode::input_error, "missing required option: --map");
    }
    if (!config.sample_sheet_path.empty() && !config.direct_samples.empty()) {
        throw Iso2GeneError(
            ExitCode::input_error,
            "use either --sample-sheet or direct sample=path inputs, not both"
        );
    }
    if (config.sample_sheet_path.empty() && config.direct_samples.empty()) {
        throw Iso2GeneError(
            ExitCode::input_error,
            "missing samples; use --sample-sheet or direct sample=path inputs"
        );
    }

    std::unordered_set<std::string> names;
    for (const SampleInput& sample : config.direct_samples) {
        if (sample.name.empty() || sample.path.empty()) {
            throw Iso2GeneError(ExitCode::input_error, "sample name and path must be non-empty");
        }
        if (!names.insert(sample.name).second) {
            throw Iso2GeneError(
                ExitCode::input_error,
                "duplicate sample name in direct inputs: " + sample.name
            );
        }
    }
}

void validate_config(const Config& config) {
    if (config.command == Command::counts || config.command == Command::txout) {
        validate_quantification_config(config);
        return;
    }
    if (config.command == Command::make_map) {
        if (config.gtf_path.empty()) {
            throw Iso2GeneError(ExitCode::input_error, "missing required option: --gtf");
        }
        if (config.map_out_path.empty()) {
            throw Iso2GeneError(ExitCode::input_error, "missing required option: --out");
        }
        if (config.transcript_id_attr.empty()) {
            throw Iso2GeneError(ExitCode::input_error, "--transcript-id-attr must be non-empty");
        }
        if (config.gene_id_attr.empty()) {
            throw Iso2GeneError(ExitCode::input_error, "--gene-id-attr must be non-empty");
        }
        return;
    }
}

} // namespace

Config parse_args(int argc, char** argv) {
    Config config;

    if (argc <= 1) {
        config.show_help = true;
        return config;
    }

    std::string first = argv[1];
    if (first == "--help" || first == "-h" || first == "help") {
        config.show_help = true;
        config.command = Command::help;
        return config;
    }
    if (first == "--version" || first == "version") {
        config.command = Command::version;
        return config;
    }
    if (first != "counts" && first != "txout" && first != "make-map") {
        throw Iso2GeneError(
            ExitCode::input_error,
            "unknown command '" + first + "'; expected counts, txout, make-map, or help"
        );
    }

    if (first == "counts") {
        config.command = Command::counts;
    } else if (first == "txout") {
        config.command = Command::txout;
    } else {
        config.command = Command::make_map;
    }

    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            config.show_help = true;
            return config;
        }

        if (config.command == Command::counts || config.command == Command::txout) {
            if (arg == "--type") {
                config.input_type = require_value(argc, argv, i, arg);
            } else if (arg == "--map") {
                config.tx2gene_path = require_value(argc, argv, i, arg);
            } else if (arg == "--sample-sheet") {
                config.sample_sheet_path = require_value(argc, argv, i, arg);
            } else if (arg == "--mode") {
                const std::string mode = require_value(argc, argv, i, arg);
                if (config.command == Command::counts) {
                    config.mode = parse_mode(mode);
                } else {
                    config.transcript_mode = parse_transcript_mode(mode);
                }
            } else if (arg == "--outdir") {
                config.outdir = require_value(argc, argv, i, arg);
            } else if (arg == "--precision") {
                config.precision = parse_precision(require_value(argc, argv, i, arg));
            } else if (arg == "--ignore-version") {
                config.id_options.ignore_version = true;
            } else if (arg == "--ignore-after-bar") {
                config.id_options.ignore_after_bar = true;
            } else if (is_option(arg)) {
                throw Iso2GeneError(ExitCode::input_error, "unknown option: " + arg);
            } else {
                config.direct_samples.push_back(parse_direct_sample(arg));
            }
        } else {
            if (arg == "--gtf") {
                config.gtf_path = require_value(argc, argv, i, arg);
            } else if (arg == "--out") {
                config.map_out_path = require_value(argc, argv, i, arg);
            } else if (arg == "--transcript-id-attr") {
                config.transcript_id_attr = require_value(argc, argv, i, arg);
            } else if (arg == "--gene-id-attr") {
                config.gene_id_attr = require_value(argc, argv, i, arg);
            } else if (is_option(arg)) {
                throw Iso2GeneError(ExitCode::input_error, "unknown option: " + arg);
            } else {
                throw Iso2GeneError(
                    ExitCode::input_error,
                    "unexpected argument for make-map: " + arg
                );
            }
        }
    }

    validate_config(config);
    return config;
}

std::string mode_name(CountMode mode) {
    switch (mode) {
    case CountMode::simple_sum:
        return "simple-sum";
    case CountMode::scaled_tpm:
        return "scaled-tpm";
    case CountMode::length_scaled_tpm:
        return "length-scaled-tpm";
    }
    throw std::logic_error("unknown count mode");
}

std::string transcript_mode_name(TranscriptCountMode mode) {
    switch (mode) {
    case TranscriptCountMode::dtu_scaled_tpm:
        return "dtu-scaled-tpm";
    }
    throw std::logic_error("unknown transcript count mode");
}

std::string help_text() {
    std::ostringstream out;
    out
        << "iso2gene: convert transcript-level quantification to gene- or transcript-level matrices\n\n"
        << "Usage:\n"
        << "  iso2gene counts --type TYPE --map tx2gene.tsv --sample-sheet samples.tsv --outdir out [options]\n"
        << "  iso2gene counts --type TYPE --map tx2gene.tsv --outdir out [options] sample=quant-file ...\n"
        << "  iso2gene txout --type TYPE --map tx2gene.tsv --sample-sheet samples.tsv --outdir out [options]\n"
        << "  iso2gene txout --type TYPE --map tx2gene.tsv --outdir out [options] sample=quant-file ...\n"
        << "  iso2gene make-map --gtf annotation.gtf[.gz] --out tx2gene.tsv [options]\n"
        << "  iso2gene --version\n\n"
        << "Required for counts and txout:\n"
        << "  --map PATH             transcript-to-gene TSV; first two columns are transcript_id and gene_id\n"
        << "  --sample-sheet PATH    TSV with sample and path columns, unless direct sample=path inputs are used\n\n"
        << "Options:\n"
        << "  --type TYPE            kallisto, salmon, or rsem (default: kallisto)\n"
        << "  --mode MODE            simple-sum, scaled-tpm, or length-scaled-tpm (default)\n"
        << "                         for txout: dtu-scaled-tpm (default)\n"
        << "  --outdir DIR           output directory (default: out)\n"
        << "  --precision N          numeric precision from 1 to 17 (default: 10)\n"
        << "  --ignore-version       strip transcript suffix after first dot\n"
        << "  --ignore-after-bar     strip transcript suffix after first bar\n"
        << "  --help                 show this help\n\n"
        << "Options for make-map:\n"
        << "  --gtf PATH             plain text or gzip-compressed GTF annotation\n"
        << "  --out PATH             output tx2gene TSV path\n"
        << "  --transcript-id-attr NAME attribute name for transcript IDs (default: transcript_id)\n"
        << "  --gene-id-attr NAME    attribute name for gene IDs (default: gene_id)\n\n"
        << "Outputs for counts:\n"
        << "  gene_counts.tsv, gene_tpm.tsv, gene_length.tsv, summary.tsv, warnings.log\n\n"
        << "Outputs for txout:\n"
        << "  transcript_counts.tsv, transcript_tpm.tsv, transcript_length.tsv, transcript_gene.tsv, summary.tsv, warnings.log\n\n"
        << "Output for make-map:\n"
        << "  tx2gene TSV at --out\n";
    return out.str();
}

std::string version_text() {
    return std::string("iso2gene ") + version_string + "\n";
}

} // namespace iso2gene
