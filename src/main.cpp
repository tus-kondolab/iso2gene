#include <exception>
#include <iostream>
#include <vector>

#include "iso2gene/cli.hpp"
#include "iso2gene/error.hpp"
#include "iso2gene/kallisto.hpp"
#include "iso2gene/logging.hpp"
#include "iso2gene/quant.hpp"
#include "iso2gene/quant_reader.hpp"
#include "iso2gene/sample_sheet.hpp"
#include "iso2gene/summarize.hpp"
#include "iso2gene/tx2gene.hpp"
#include "iso2gene/write_matrix.hpp"

namespace {

int run_counts(const iso2gene::Config& config) {
    iso2gene::Logger logger;

    std::vector<iso2gene::SampleInput> samples;
    if (!config.sample_sheet_path.empty()) {
        samples = iso2gene::read_sample_sheet(config.sample_sheet_path);
    } else {
        samples = config.direct_samples;
    }

    logger.info("reading tx2gene map");
    const iso2gene::Tx2GeneMap tx2gene =
        iso2gene::read_tx2gene(config.tx2gene_path, config.id_options, logger);

    std::vector<std::vector<iso2gene::QuantRecord>> quantifications;
    quantifications.reserve(samples.size());
    for (const iso2gene::SampleInput& sample : samples) {
        logger.info("reading sample '" + sample.name + "'");
        quantifications.push_back(iso2gene::read_quantification_file(
            config.input_type,
            sample.path,
            config.id_options
        ));
    }

    logger.info("summarizing transcripts to genes");
    const iso2gene::GeneMatrices matrices =
        iso2gene::summarize_to_gene(samples, quantifications, tx2gene, config.mode, logger);

    logger.info("writing outputs");
    iso2gene::write_outputs(config, matrices, logger);
    return static_cast<int>(iso2gene::ExitCode::success);
}

} // namespace

int main(int argc, char** argv) {
    try {
        const iso2gene::Config config = iso2gene::parse_args(argc, argv);
        if (config.show_help || config.command == iso2gene::Command::help) {
            std::cout << iso2gene::help_text();
            return static_cast<int>(iso2gene::ExitCode::success);
        }

        switch (config.command) {
        case iso2gene::Command::counts:
            return run_counts(config);
        case iso2gene::Command::help:
            std::cout << iso2gene::help_text();
            return static_cast<int>(iso2gene::ExitCode::success);
        }
    } catch (const iso2gene::Iso2GeneError& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return static_cast<int>(error.code());
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return static_cast<int>(iso2gene::ExitCode::internal_error);
    } catch (...) {
        std::cerr << "ERROR: unknown internal error\n";
        return static_cast<int>(iso2gene::ExitCode::internal_error);
    }

    std::cerr << "ERROR: unreachable command state\n";
    return static_cast<int>(iso2gene::ExitCode::internal_error);
}
