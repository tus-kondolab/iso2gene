#include <cmath>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include "iso2gene/id.hpp"
#include "iso2gene/kallisto.hpp"
#include "iso2gene/logging.hpp"
#include "iso2gene/quant.hpp"
#include "iso2gene/quant_reader.hpp"
#include "iso2gene/rsem.hpp"
#include "iso2gene/salmon.hpp"
#include "iso2gene/sample_sheet.hpp"
#include "iso2gene/summarize.hpp"
#include "iso2gene/tsv.hpp"
#include "iso2gene/tx2gene.hpp"
#include "iso2gene/error.hpp"

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_near(double actual, double expected, const std::string& message) {
    const double diff = std::fabs(actual - expected);
    const double scale = std::max(1.0, std::fabs(expected));
    if (diff > 1e-6 * scale) {
        throw std::runtime_error(
            message + ": expected " + std::to_string(expected)
            + ", got " + std::to_string(actual)
        );
    }
}

std::size_t gene_index(const iso2gene::GeneMatrices& matrices, const std::string& gene_id) {
    for (std::size_t i = 0; i < matrices.gene_ids.size(); ++i) {
        if (matrices.gene_ids[i] == gene_id) {
            return i;
        }
    }
    throw std::runtime_error("missing gene in output: " + gene_id);
}

void test_split_tsv_line() {
    const std::vector<std::string> fields = iso2gene::split_tsv_line("a\t\tb\t");
    require(fields.size() == 4, "split_tsv_line preserves trailing empty fields");
    require(fields[0] == "a", "split field 0");
    require(fields[1].empty(), "split field 1 empty");
    require(fields[2] == "b", "split field 2");
    require(fields[3].empty(), "split field 3 empty");
}

void test_id_normalization() {
    iso2gene::IdOptions options;
    options.ignore_after_bar = true;
    options.ignore_version = true;
    require(
        iso2gene::normalize_transcript_id("ENST0001.5|ENSG0001|x", options) == "ENST0001",
        "id normalization applies bar then version stripping"
    );
}

void test_readers() {
    iso2gene::Logger logger;
    const iso2gene::IdOptions id_options;

    const iso2gene::Tx2GeneMap map =
        iso2gene::read_tx2gene("tests/data/tx2gene.tsv", id_options, logger);
    require(map.tx_to_gene.size() == 4, "tx2gene size");
    require(map.gene_ids.size() == 3, "gene id count");
    require(map.tx_to_gene.at("tx2") == "geneA", "tx2gene lookup");

    const std::vector<iso2gene::QuantRecord> records =
        iso2gene::read_kallisto("tests/data/sample1_abundance.tsv", id_options);
    require(records.size() == 3, "kallisto record count");
    require(records[0].transcript_id == "tx1", "kallisto transcript id");
    require_near(records[2].tpm, 200000.0, "kallisto tpm");

    const std::vector<iso2gene::QuantRecord> salmon_records =
        iso2gene::read_salmon("tests/data/salmon_sample1_quant.sf", id_options);
    require(salmon_records.size() == 3, "salmon record count");
    require(salmon_records[0].transcript_id == "tx1", "salmon transcript id");
    require_near(salmon_records[0].tpm, 600000.0, "salmon scientific notation tpm");
    require_near(salmon_records[2].est_counts, 10.0, "salmon scientific notation NumReads");

    const std::vector<iso2gene::QuantRecord> rsem_records =
        iso2gene::read_rsem("tests/data/rsem_sample1_isoforms.results", id_options);
    require(rsem_records.size() == 3, "rsem record count");
    require(rsem_records[0].transcript_id == "tx1", "rsem transcript id");
    require_near(rsem_records[2].eff_length, 900.0, "rsem effective length");

    const std::vector<iso2gene::QuantRecord> rsem_zero_records =
        iso2gene::read_rsem("tests/data/rsem_zero_eff_isoforms.results", id_options);
    require_near(rsem_zero_records[0].eff_length, 1.0, "rsem zero effective length is clamped like tximport");

    const std::vector<iso2gene::SampleInput> samples =
        iso2gene::read_sample_sheet("tests/data/sample_sheet.tsv");
    require(samples.size() == 2, "sample sheet count");
    require(samples[1].name == "s2", "sample sheet sample name");
}

void test_quantification_dispatch() {
    const iso2gene::IdOptions id_options;

    const std::vector<iso2gene::QuantRecord> kallisto_records =
        iso2gene::read_quantification_file(
            "kallisto",
            "tests/data/sample1_abundance.tsv",
            id_options
        );
    require_near(kallisto_records[0].est_counts, 30.0, "dispatch kallisto reader");

    const std::vector<iso2gene::QuantRecord> salmon_records =
        iso2gene::read_quantification_file(
            "salmon",
            "tests/data/salmon_sample1_quant.sf",
            id_options
        );
    require_near(salmon_records[0].tpm, 600000.0, "dispatch salmon reader");

    const std::vector<iso2gene::QuantRecord> rsem_records =
        iso2gene::read_quantification_file(
            "rsem",
            "tests/data/rsem_sample1_isoforms.results",
            id_options
        );
    require_near(rsem_records[0].eff_length, 90.0, "dispatch rsem reader");

    bool saw_error = false;
    try {
        (void)iso2gene::read_quantification_file(
            "unknown",
            "tests/data/sample1_abundance.tsv",
            id_options
        );
    } catch (const iso2gene::Iso2GeneError& error) {
        saw_error = error.code() == iso2gene::ExitCode::input_error;
    }
    require(saw_error, "dispatch rejects unknown input type");
}

iso2gene::GeneMatrices summarize_fixture(iso2gene::CountMode mode) {
    iso2gene::Logger logger;
    const iso2gene::IdOptions id_options;
    const iso2gene::Tx2GeneMap map =
        iso2gene::read_tx2gene("tests/data/tx2gene.tsv", id_options, logger);
    const std::vector<iso2gene::SampleInput> samples{
        {"s1", "tests/data/sample1_abundance.tsv"},
        {"s2", "tests/data/sample2_abundance.tsv"}
    };
    std::vector<std::vector<iso2gene::QuantRecord>> quantifications;
    quantifications.push_back(iso2gene::read_kallisto(samples[0].path, id_options));
    quantifications.push_back(iso2gene::read_kallisto(samples[1].path, id_options));
    return iso2gene::summarize_to_gene(samples, quantifications, map, mode, logger);
}

using QuantReader = std::vector<iso2gene::QuantRecord> (*)(
    const std::string&,
    const iso2gene::IdOptions&
);

iso2gene::GeneMatrices summarize_fixture_with_reader(
    QuantReader reader,
    const std::string& sample1_path,
    const std::string& sample2_path,
    iso2gene::CountMode mode
) {
    iso2gene::Logger logger;
    const iso2gene::IdOptions id_options;
    const iso2gene::Tx2GeneMap map =
        iso2gene::read_tx2gene("tests/data/tx2gene.tsv", id_options, logger);
    const std::vector<iso2gene::SampleInput> samples{
        {"s1", sample1_path},
        {"s2", sample2_path}
    };
    std::vector<std::vector<iso2gene::QuantRecord>> quantifications;
    quantifications.push_back(reader(samples[0].path, id_options));
    quantifications.push_back(reader(samples[1].path, id_options));
    return iso2gene::summarize_to_gene(samples, quantifications, map, mode, logger);
}

void test_simple_sum() {
    const iso2gene::GeneMatrices matrices =
        summarize_fixture(iso2gene::CountMode::simple_sum);
    require(matrices.gene_ids.size() == 2, "unobserved tx2gene genes are not emitted");
    require(matrices.gene_ids[0] == "geneA", "gene order geneA");
    require(matrices.gene_ids[1] == "geneB", "gene order geneB");

    const std::size_t geneA = gene_index(matrices, "geneA");
    const std::size_t geneB = gene_index(matrices, "geneB");

    require_near(matrices.counts(geneA, 0), 40.0, "simple geneA s1 counts");
    require_near(matrices.counts(geneA, 1), 40.0, "simple geneA s2 counts");
    require_near(matrices.counts(geneB, 0), 10.0, "simple geneB s1 counts");
    require_near(matrices.counts(geneB, 1), 60.0, "simple geneB s2 counts");

    require_near(matrices.tpm(geneA, 0), 800000.0, "geneA s1 tpm");
    require_near(matrices.tpm(geneA, 1), 500000.0, "geneA s2 tpm");
    require_near(matrices.length(geneA, 0), 112.5, "geneA s1 length");
    require_near(matrices.length(geneA, 1), 160.0, "geneA s2 length");
    require_near(matrices.length(geneB, 0), 900.0, "geneB s1 length");
}

void test_scaled_tpm() {
    const iso2gene::GeneMatrices matrices =
        summarize_fixture(iso2gene::CountMode::scaled_tpm);
    const std::size_t geneA = gene_index(matrices, "geneA");
    const std::size_t geneB = gene_index(matrices, "geneB");
    require_near(matrices.counts(geneA, 0), 40.0, "scaled geneA s1 counts");
    require_near(matrices.counts(geneB, 0), 10.0, "scaled geneB s1 counts");
    require_near(matrices.counts(geneA, 1), 50.0, "scaled geneA s2 counts");
    require_near(matrices.counts(geneB, 1), 50.0, "scaled geneB s2 counts");
}

void test_length_scaled_tpm() {
    const iso2gene::GeneMatrices matrices =
        summarize_fixture(iso2gene::CountMode::length_scaled_tpm);
    const std::size_t geneA = gene_index(matrices, "geneA");
    const std::size_t geneB = gene_index(matrices, "geneB");
    require_near(matrices.counts(geneA, 0), 18.2274247492, "length-scaled geneA s1 counts");
    require_near(matrices.counts(geneB, 0), 31.7725752508, "length-scaled geneB s1 counts");
    require_near(matrices.counts(geneA, 1), 12.5431530495, "length-scaled geneA s2 counts");
    require_near(matrices.counts(geneB, 1), 87.4568469505, "length-scaled geneB s2 counts");
}

void assert_fixture_matrix_values(const iso2gene::GeneMatrices& matrices, const std::string& label) {
    const std::size_t geneA = gene_index(matrices, "geneA");
    const std::size_t geneB = gene_index(matrices, "geneB");
    require_near(matrices.counts(geneA, 0), 18.2274247492, label + " geneA s1 counts");
    require_near(matrices.counts(geneB, 0), 31.7725752508, label + " geneB s1 counts");
    require_near(matrices.counts(geneA, 1), 12.5431530495, label + " geneA s2 counts");
    require_near(matrices.counts(geneB, 1), 87.4568469505, label + " geneB s2 counts");
    require_near(matrices.tpm(geneA, 0), 800000.0, label + " geneA s1 tpm");
    require_near(matrices.length(geneA, 0), 112.5, label + " geneA s1 length");
}

void test_reader_formats_share_summarization() {
    const iso2gene::GeneMatrices salmon =
        summarize_fixture_with_reader(
            iso2gene::read_salmon,
            "tests/data/salmon_sample1_quant.sf",
            "tests/data/salmon_sample2_quant.sf",
            iso2gene::CountMode::length_scaled_tpm
        );
    assert_fixture_matrix_values(salmon, "salmon");

    const iso2gene::GeneMatrices rsem =
        summarize_fixture_with_reader(
            iso2gene::read_rsem,
            "tests/data/rsem_sample1_isoforms.results",
            "tests/data/rsem_sample2_isoforms.results",
            iso2gene::CountMode::length_scaled_tpm
        );
    assert_fixture_matrix_values(rsem, "rsem");
}

void test_scaled_tpm_uses_mapped_denominators() {
    iso2gene::Logger logger;
    const iso2gene::IdOptions id_options;
    const iso2gene::Tx2GeneMap map =
        iso2gene::read_tx2gene("tests/data/tx2gene_unmapped.tsv", id_options, logger);
    const std::vector<iso2gene::SampleInput> samples{
        {"s1", "tests/data/sample_unmapped_abundance.tsv"}
    };
    std::vector<std::vector<iso2gene::QuantRecord>> quantifications;
    quantifications.push_back(iso2gene::read_kallisto(samples[0].path, id_options));

    const iso2gene::GeneMatrices scaled =
        iso2gene::summarize_to_gene(samples, quantifications, map, iso2gene::CountMode::scaled_tpm, logger);
    require_near(scaled.counts(gene_index(scaled, "geneA"), 0), 10.0, "mapped scaledTPM geneA");
    require_near(scaled.counts(gene_index(scaled, "geneB"), 0), 30.0, "mapped scaledTPM geneB");

    const iso2gene::GeneMatrices length_scaled =
        iso2gene::summarize_to_gene(samples, quantifications, map, iso2gene::CountMode::length_scaled_tpm, logger);
    require_near(length_scaled.counts(gene_index(length_scaled, "geneA"), 0), 5.7142857143, "mapped lengthScaledTPM geneA");
    require_near(length_scaled.counts(gene_index(length_scaled, "geneB"), 0), 34.2857142857, "mapped lengthScaledTPM geneB");
}

void test_zero_abundance_length_replacement() {
    iso2gene::Logger logger;
    const iso2gene::IdOptions id_options;
    const iso2gene::Tx2GeneMap map =
        iso2gene::read_tx2gene("tests/data/tx2gene_zero.tsv", id_options, logger);
    const std::vector<iso2gene::SampleInput> samples{
        {"s1", "tests/data/zero_s1_abundance.tsv"},
        {"s2", "tests/data/zero_s2_abundance.tsv"}
    };
    std::vector<std::vector<iso2gene::QuantRecord>> quantifications;
    quantifications.push_back(iso2gene::read_kallisto(samples[0].path, id_options));
    quantifications.push_back(iso2gene::read_kallisto(samples[1].path, id_options));

    const iso2gene::GeneMatrices matrices =
        iso2gene::summarize_to_gene(samples, quantifications, map, iso2gene::CountMode::simple_sum, logger);
    require_near(matrices.length(gene_index(matrices, "geneA"), 0), 185.0, "geneA missing length replaced from other sample");
    require_near(matrices.length(gene_index(matrices, "geneA"), 1), 185.0, "geneA observed weighted length");
    require_near(matrices.length(gene_index(matrices, "geneC"), 0), 600.0, "geneC all-zero length replacement s1");
    require_near(matrices.length(gene_index(matrices, "geneC"), 1), 600.0, "geneC all-zero length replacement s2");
}

void test_rsem_zero_effective_length_summarization() {
    iso2gene::Logger logger;
    const iso2gene::IdOptions id_options;
    const iso2gene::Tx2GeneMap map =
        iso2gene::read_tx2gene("tests/data/tx2gene_rsem_zero_eff.tsv", id_options, logger);
    const std::vector<iso2gene::SampleInput> samples{
        {"s1", "tests/data/rsem_zero_eff_isoforms.results"}
    };
    std::vector<std::vector<iso2gene::QuantRecord>> quantifications;
    quantifications.push_back(iso2gene::read_rsem(samples[0].path, id_options));

    const iso2gene::GeneMatrices matrices =
        iso2gene::summarize_to_gene(samples, quantifications, map, iso2gene::CountMode::simple_sum, logger);
    const std::size_t geneZ = gene_index(matrices, "geneZ");
    require_near(matrices.counts(geneZ, 0), 0.0, "rsem zero eff length count");
    require_near(matrices.tpm(geneZ, 0), 0.0, "rsem zero eff length tpm");
    require_near(matrices.length(geneZ, 0), 100.5, "rsem all-zero gene length replacement clamps zero effective length");
}

} // namespace

int main() {
    try {
        test_split_tsv_line();
        test_id_normalization();
        test_readers();
        test_quantification_dispatch();
        test_simple_sum();
        test_scaled_tpm();
        test_length_scaled_tpm();
        test_reader_formats_share_summarization();
        test_scaled_tpm_uses_mapped_denominators();
        test_zero_abundance_length_replacement();
        test_rsem_zero_effective_length_summarization();
    } catch (const std::exception& error) {
        std::cerr << "TEST FAILED: " << error.what() << '\n';
        return 1;
    }

    std::cout << "all tests passed\n";
    return 0;
}
