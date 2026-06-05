#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "iso2gene/cli.hpp"
#include "iso2gene/gtf.hpp"
#include "iso2gene/id.hpp"
#include "iso2gene/kallisto.hpp"
#include "iso2gene/logging.hpp"
#include "iso2gene/quant.hpp"
#include "iso2gene/quant_reader.hpp"
#include "iso2gene/rsem.hpp"
#include "iso2gene/salmon.hpp"
#include "iso2gene/sample_sheet.hpp"
#include "iso2gene/summarize.hpp"
#include "iso2gene/text_reader.hpp"
#include "iso2gene/tsv.hpp"
#include "iso2gene/tx2gene.hpp"
#include "iso2gene/write_matrix.hpp"
#include "iso2gene/error.hpp"
#include "iso2gene/version.hpp"

#include "miniz.h"

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

std::string temp_output_path(const std::string& filename) {
    return (std::filesystem::temp_directory_path() / filename).string();
}

std::string read_binary_file(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("failed to open fixture: " + path);
    }
    return std::string(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>()
    );
}

void write_binary_file(const std::string& path, const std::string& data) {
    std::ofstream output(path, std::ios::binary);
    if (!output) {
        throw std::runtime_error("failed to open temp output: " + path);
    }
    output.write(data.data(), static_cast<std::streamsize>(data.size()));
    if (!output) {
        throw std::runtime_error("failed to write temp output: " + path);
    }
}

void append_u32_le(std::string& data, const std::uint32_t value) {
    data.push_back(static_cast<char>(value & 0xFFU));
    data.push_back(static_cast<char>((value >> 8U) & 0xFFU));
    data.push_back(static_cast<char>((value >> 16U) & 0xFFU));
    data.push_back(static_cast<char>((value >> 24U) & 0xFFU));
}

struct GzipTestHeaderOptions {
    bool extra = false;
    bool name = false;
    bool comment = false;
    bool header_crc = false;
};

std::string gzip_bytes(
    const std::string& data,
    const bool corrupt_crc = false,
    const bool truncate = false,
    const GzipTestHeaderOptions& header_options = {}
) {
    std::size_t compressed_size = 0;
    void* compressed = tdefl_compress_mem_to_heap(
        data.data(),
        data.size(),
        &compressed_size,
        TDEFL_DEFAULT_MAX_PROBES
    );
    if (compressed == nullptr) {
        throw std::runtime_error("failed to gzip test fixture");
    }

    std::string output;
    output.push_back(static_cast<char>(0x1F));
    output.push_back(static_cast<char>(0x8B));
    output.push_back(static_cast<char>(0x08));
    unsigned char flags = 0;
    if (header_options.extra) {
        flags |= 0x04U;
    }
    if (header_options.name) {
        flags |= 0x08U;
    }
    if (header_options.comment) {
        flags |= 0x10U;
    }
    if (header_options.header_crc) {
        flags |= 0x02U;
    }
    output.push_back(static_cast<char>(flags));
    output.append(4, '\0');
    output.push_back(static_cast<char>(0x00));
    output.push_back(static_cast<char>(0xFF));
    if (header_options.extra) {
        output.push_back(static_cast<char>(0x04));
        output.push_back(static_cast<char>(0x00));
        output.append("I2GZ", 4);
    }
    if (header_options.name) {
        output.append("sample.tsv");
        output.push_back('\0');
    }
    if (header_options.comment) {
        output.append("iso2gene gzip optional header test");
        output.push_back('\0');
    }
    if (header_options.header_crc) {
        output.push_back('\0');
        output.push_back('\0');
    }
    output.append(
        static_cast<const char*>(compressed),
        static_cast<const char*>(compressed) + compressed_size
    );
    mz_free(compressed);

    std::uint32_t crc = static_cast<std::uint32_t>(
        mz_crc32(static_cast<mz_ulong>(MZ_CRC32_INIT),
                 reinterpret_cast<const unsigned char*>(data.data()),
                 data.size())
    );
    if (corrupt_crc) {
        crc ^= 0xFFFFFFFFU;
    }
    append_u32_le(output, crc);
    append_u32_le(output, static_cast<std::uint32_t>(data.size()));

    if (truncate && output.size() > 3) {
        output.resize(output.size() - 3);
    }
    return output;
}

std::string gzip_file_to_temp(const std::string& input_path, const std::string& filename) {
    const std::string output_path = temp_output_path(filename);
    write_binary_file(output_path, gzip_bytes(read_binary_file(input_path)));
    return output_path;
}

std::vector<std::string> read_all_lines(const std::string& path) {
    iso2gene::TextReader reader(path);
    std::vector<std::string> lines;
    std::string line;
    while (reader.read_line(line)) {
        lines.push_back(line);
    }
    return lines;
}

void drain_text_reader(const std::string& path) {
    iso2gene::TextReader reader(path);
    std::string line;
    while (reader.read_line(line)) {
    }
}

void require_lf_only_file(const std::string& path, const std::string& label) {
    const std::string contents = read_binary_file(path);
    require(contents.find('\n') != std::string::npos, label + " contains LF");
    require(contents.find('\r') == std::string::npos, label + " does not contain CR");
}

void test_make_map_cli_parse() {
    char arg0[] = "iso2gene";
    char arg1[] = "make-map";
    char arg2[] = "--gtf";
    char arg3[] = "annotation.gtf";
    char arg4[] = "--out";
    char arg5[] = "tx2gene.tsv";
    char arg6[] = "--transcript-id-attr";
    char arg7[] = "transcript";
    char arg8[] = "--gene-id-attr";
    char arg9[] = "gene";
    char* argv[] = {arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9};

    const iso2gene::Config config = iso2gene::parse_args(10, argv);
    require(config.command == iso2gene::Command::make_map, "make-map command parse");
    require(config.gtf_path == "annotation.gtf", "make-map --gtf parse");
    require(config.map_out_path == "tx2gene.tsv", "make-map --out parse");
    require(config.transcript_id_attr == "transcript", "make-map transcript attr parse");
    require(config.gene_id_attr == "gene", "make-map gene attr parse");
}

void test_version_cli_parse() {
    char arg0[] = "iso2gene";
    char arg1[] = "--version";
    char* argv[] = {arg0, arg1};

    const iso2gene::Config config = iso2gene::parse_args(2, argv);
    require(config.command == iso2gene::Command::version, "version command parse");
    require(std::string(iso2gene::version_string) == "1.0.0", "version string");
    require(iso2gene::version_text() == "iso2gene 1.0.0\n", "version text");
}

void test_make_map_basic_gtf() {
    iso2gene::Logger logger;
    const std::string out_path = temp_output_path("iso2gene_annotation_basic_tx2gene.tsv");
    const iso2gene::GtfMapOptions options{
        "tests/data/annotation_basic.gtf",
        out_path,
        "transcript_id",
        "gene_id"
    };

    const iso2gene::GtfMapStats stats = iso2gene::make_tx2gene_from_gtf(options, logger);
    require(stats.mappings_written == 3, "GTF basic mappings written");
    require(stats.duplicate_same_gene_rows == 1, "GTF duplicate same-gene rows");
    require(stats.missing_transcript_id_rows == 1, "GTF gene rows without transcript are skipped");

    const iso2gene::IdOptions id_options;
    const iso2gene::Tx2GeneMap map = iso2gene::read_tx2gene(out_path, id_options, logger);
    require(map.tx_to_gene.size() == 3, "generated tx2gene size");
    require(map.tx_to_gene.at("tx1") == "geneA", "generated tx1 mapping");
    require(map.tx_to_gene.at("tx2") == "geneB", "generated tx2 mapping");
    require(map.tx_to_gene.at("tx3") == "geneC", "generated tx3 mapping");
}

void test_make_map_gencode_attribute_edges() {
    iso2gene::Logger logger;
    const std::string out_path = temp_output_path("iso2gene_annotation_gencode_tx2gene.tsv");
    const iso2gene::GtfMapOptions options{
        "tests/data/annotation_gencode_fragment.gtf",
        out_path,
        "transcript_id",
        "gene_id"
    };

    const iso2gene::GtfMapStats stats = iso2gene::make_tx2gene_from_gtf(options, logger);
    require(stats.mappings_written == 1, "GENCODE-like fixture mapping count");
    require(stats.empty_transcript_id_rows == 1, "empty transcript_id is skipped");
    require(stats.missing_gene_id_rows == 1, "similar gene keys are not gene_id");

    const iso2gene::IdOptions id_options;
    const iso2gene::Tx2GeneMap map = iso2gene::read_tx2gene(out_path, id_options, logger);
    require(map.tx_to_gene.size() == 1, "GENCODE-like output size");
    require(
        map.tx_to_gene.at("ENST000001.1") == "ENSG000001.1",
        "GENCODE-like exact gene_id mapping"
    );
    require(
        map.tx_to_gene.find("ENSTREF.1") == map.tx_to_gene.end(),
        "ref_gene_id is not mistaken for gene_id"
    );
}

void test_make_map_rejects_conflicts_and_bad_attributes() {
    iso2gene::Logger logger;

    bool saw_conflict = false;
    try {
        const iso2gene::GtfMapOptions options{
            "tests/data/annotation_conflict.gtf",
            temp_output_path("iso2gene_annotation_conflict_tx2gene.tsv"),
            "transcript_id",
            "gene_id"
        };
        (void)iso2gene::make_tx2gene_from_gtf(options, logger);
    } catch (const iso2gene::Iso2GeneError& error) {
        saw_conflict = error.code() == iso2gene::ExitCode::parse_error
            && std::string(error.what()).find("maps to multiple genes") != std::string::npos;
    }
    require(saw_conflict, "GTF conflict mapping is rejected");

    bool saw_duplicate_attr = false;
    try {
        const iso2gene::GtfMapOptions options{
            "tests/data/annotation_bad_duplicate_attr.gtf",
            temp_output_path("iso2gene_annotation_bad_duplicate_attr_tx2gene.tsv"),
            "transcript_id",
            "gene_id"
        };
        (void)iso2gene::make_tx2gene_from_gtf(options, logger);
    } catch (const iso2gene::Iso2GeneError& error) {
        saw_duplicate_attr = error.code() == iso2gene::ExitCode::parse_error
            && std::string(error.what()).find("appears multiple times") != std::string::npos;
    }
    require(saw_duplicate_attr, "conflicting duplicate transcript_id attribute is rejected");

    bool saw_gff3_hint = false;
    try {
        const iso2gene::GtfMapOptions options{
            "tests/data/annotation_gff3.gtf",
            temp_output_path("iso2gene_annotation_gff3_tx2gene.tsv"),
            "transcript_id",
            "gene_id"
        };
        (void)iso2gene::make_tx2gene_from_gtf(options, logger);
    } catch (const iso2gene::Iso2GeneError& error) {
        saw_gff3_hint = error.code() == iso2gene::ExitCode::parse_error
            && std::string(error.what()).find("GFF3 is not supported") != std::string::npos;
    }
    require(saw_gff3_hint, "GFF3-like attributes get a helpful error");
}

void test_text_reader_line_endings_and_bom() {
    const std::vector<std::string> variants{
        "a\tb\nc\td\n",
        "a\tb\r\nc\td\r\n",
        "a\tb\rc\td\r",
        "a\tb\nc\td"
    };

    for (std::size_t i = 0; i < variants.size(); ++i) {
        const std::string plain_path =
            temp_output_path("iso2gene_line_endings_" + std::to_string(i) + ".txt");
        write_binary_file(plain_path, variants[i]);
        const std::vector<std::string> plain_lines = read_all_lines(plain_path);
        require(plain_lines.size() == 2, "plain line ending logical line count");
        require(plain_lines[0] == "a\tb", "plain line ending first line");
        require(plain_lines[1] == "c\td", "plain line ending second line");

        const std::string gzip_path = plain_path + ".gz";
        write_binary_file(gzip_path, gzip_bytes(variants[i]));
        const std::vector<std::string> gzip_lines = read_all_lines(gzip_path);
        require(gzip_lines == plain_lines, "gzip line endings match plain lines");
    }

    const std::string bom_sheet =
        std::string("\xEF\xBB\xBF", 3) + "sample\tpath\ns1\tquant.sf\n";
    const std::string plain_bom = temp_output_path("iso2gene_bom_sample_sheet.tsv");
    write_binary_file(plain_bom, bom_sheet);
    const std::vector<iso2gene::SampleInput> plain_samples =
        iso2gene::read_sample_sheet(plain_bom);
    require(plain_samples.size() == 1, "plain BOM sample sheet count");
    require(plain_samples[0].name == "s1", "plain BOM sample sheet name");

    const std::string gzip_bom = plain_bom + ".gz";
    write_binary_file(gzip_bom, gzip_bytes(bom_sheet));
    const std::vector<iso2gene::SampleInput> gzip_samples =
        iso2gene::read_sample_sheet(gzip_bom);
    require(gzip_samples.size() == 1, "gzip BOM sample sheet count");
    require(gzip_samples[0].path == "quant.sf", "gzip BOM sample sheet path");
}

void test_external_gzip_fixture() {
    // Generated with external gzip, see tests/data/gzip/README.md.
    iso2gene::Logger logger;
    const iso2gene::IdOptions id_options;
    const iso2gene::Tx2GeneMap map = iso2gene::read_tx2gene(
        "tests/data/gzip/tx2gene_external.tsv.gz",
        id_options,
        logger
    );
    require(map.tx_to_gene.size() == 2, "external gzip tx2gene size");
    require(map.tx_to_gene.at("tx_ext1") == "geneExtA", "external gzip tx_ext1");
    require(map.tx_to_gene.at("tx_ext2") == "geneExtB", "external gzip tx_ext2");
}

void test_gzip_optional_header_fields() {
    GzipTestHeaderOptions options;
    options.extra = true;
    options.name = true;
    options.comment = true;
    options.header_crc = true;

    const std::string path = temp_output_path("iso2gene_optional_header.tsv.gz");
    write_binary_file(
        path,
        gzip_bytes("transcript_id\tgene_id\nopt_tx1\topt_gene1\n", false, false, options)
    );

    iso2gene::Logger logger;
    const iso2gene::IdOptions id_options;
    const iso2gene::Tx2GeneMap map = iso2gene::read_tx2gene(path, id_options, logger);
    require(map.tx_to_gene.size() == 1, "gzip optional header tx2gene size");
    require(map.tx_to_gene.at("opt_tx1") == "opt_gene1", "gzip optional header tx2gene mapping");
}

void test_gzip_inputs_match_plain() {
    iso2gene::Logger logger;
    const iso2gene::IdOptions id_options;

    const std::string kallisto_s1 =
        gzip_file_to_temp("tests/data/sample1_abundance.tsv", "iso2gene_sample1_abundance.tsv.gz");
    const std::string kallisto_s2 =
        gzip_file_to_temp("tests/data/sample2_abundance.tsv", "iso2gene_sample2_abundance.tsv.gz");
    const iso2gene::GeneMatrices kallisto =
        summarize_fixture_with_reader(
            iso2gene::read_kallisto,
            kallisto_s1,
            kallisto_s2,
            iso2gene::CountMode::length_scaled_tpm
        );
    assert_fixture_matrix_values(kallisto, "kallisto gzip");

    const std::string salmon_s1 =
        gzip_file_to_temp("tests/data/salmon_sample1_quant.sf", "iso2gene_salmon_sample1_quant.sf.gz");
    const std::string salmon_s2 =
        gzip_file_to_temp("tests/data/salmon_sample2_quant.sf", "iso2gene_salmon_sample2_quant.sf.gz");
    const iso2gene::GeneMatrices salmon =
        summarize_fixture_with_reader(
            iso2gene::read_salmon,
            salmon_s1,
            salmon_s2,
            iso2gene::CountMode::length_scaled_tpm
        );
    assert_fixture_matrix_values(salmon, "salmon gzip");

    const std::string rsem_s1 = gzip_file_to_temp(
        "tests/data/rsem_sample1_isoforms.results",
        "iso2gene_rsem_sample1_isoforms.results.gz"
    );
    const std::string rsem_s2 = gzip_file_to_temp(
        "tests/data/rsem_sample2_isoforms.results",
        "iso2gene_rsem_sample2_isoforms.results.gz"
    );
    const iso2gene::GeneMatrices rsem =
        summarize_fixture_with_reader(
            iso2gene::read_rsem,
            rsem_s1,
            rsem_s2,
            iso2gene::CountMode::length_scaled_tpm
        );
    assert_fixture_matrix_values(rsem, "rsem gzip");

    const std::string tx2gene_gz =
        gzip_file_to_temp("tests/data/tx2gene.tsv", "iso2gene_tx2gene.tsv.gz");
    const iso2gene::Tx2GeneMap map = iso2gene::read_tx2gene(tx2gene_gz, id_options, logger);
    require(map.tx_to_gene.size() == 4, "gzip tx2gene size");

    const std::string sample_sheet_gz =
        gzip_file_to_temp("tests/data/sample_sheet.tsv", "iso2gene_sample_sheet.tsv.gz");
    const std::vector<iso2gene::SampleInput> samples =
        iso2gene::read_sample_sheet(sample_sheet_gz);
    require(samples.size() == 2, "gzip sample sheet count");

    const std::string gtf_gz =
        gzip_file_to_temp("tests/data/annotation_basic.gtf", "iso2gene_annotation_basic.gtf.gz");
    const std::string out_path = temp_output_path("iso2gene_annotation_basic_gzip_tx2gene.tsv");
    const iso2gene::GtfMapOptions options{gtf_gz, out_path, "transcript_id", "gene_id"};
    const iso2gene::GtfMapStats stats = iso2gene::make_tx2gene_from_gtf(options, logger);
    require(stats.mappings_written == 3, "gzip GTF mappings written");
    const iso2gene::Tx2GeneMap generated = iso2gene::read_tx2gene(out_path, id_options, logger);
    require(generated.tx_to_gene.at("tx3") == "geneC", "gzip GTF generated tx3 mapping");
}

void test_gzip_concatenated_members_and_errors() {
    const std::string concatenated = temp_output_path("iso2gene_concatenated.txt.gz");
    write_binary_file(concatenated, gzip_bytes("a\n") + gzip_bytes("b\n"));
    const std::vector<std::string> lines = read_all_lines(concatenated);
    require(lines.size() == 2, "concatenated gzip line count");
    require(lines[0] == "a" && lines[1] == "b", "concatenated gzip lines");

    const std::string not_gzip = temp_output_path("iso2gene_not_gzip.txt.gz");
    write_binary_file(not_gzip, "not gzip\n");
    bool saw_invalid_header = false;
    try {
        drain_text_reader(not_gzip);
    } catch (const iso2gene::Iso2GeneError& error) {
        saw_invalid_header = error.code() == iso2gene::ExitCode::io_error
            && std::string(error.what()).find("invalid gzip header") != std::string::npos;
    }
    require(saw_invalid_header, "invalid gzip header is rejected");

    const std::string truncated = temp_output_path("iso2gene_truncated.txt.gz");
    write_binary_file(truncated, gzip_bytes("a\n", false, true));
    bool saw_truncated = false;
    try {
        drain_text_reader(truncated);
    } catch (const iso2gene::Iso2GeneError& error) {
        saw_truncated = error.code() == iso2gene::ExitCode::io_error
            && std::string(error.what()).find("truncated") != std::string::npos;
    }
    require(saw_truncated, "truncated gzip is rejected");

    const std::string crc_mismatch = temp_output_path("iso2gene_crc_mismatch.txt.gz");
    write_binary_file(crc_mismatch, gzip_bytes("a\n", true, false));
    bool saw_crc_mismatch = false;
    try {
        drain_text_reader(crc_mismatch);
    } catch (const iso2gene::Iso2GeneError& error) {
        saw_crc_mismatch = error.code() == iso2gene::ExitCode::io_error
            && std::string(error.what()).find("CRC mismatch") != std::string::npos;
    }
    require(saw_crc_mismatch, "gzip CRC mismatch is rejected");
}

void test_file_outputs_use_lf() {
    const std::string outdir = temp_output_path("iso2gene_lf_outputs");

    iso2gene::Config config;
    config.outdir = outdir;
    config.input_type = "kallisto";
    config.mode = iso2gene::CountMode::length_scaled_tpm;
    config.precision = 17;

    iso2gene::GeneMatrices matrices;
    matrices.gene_ids = {"gene_lf"};
    matrices.sample_names = {"sample_lf"};
    matrices.counts = iso2gene::Matrix<double>(1, 1);
    matrices.tpm = iso2gene::Matrix<double>(1, 1);
    matrices.length = iso2gene::Matrix<double>(1, 1);
    matrices.counts(0, 0) = 1.5;
    matrices.tpm(0, 0) = 200000.0;
    matrices.length(0, 0) = 350.0;
    matrices.total_records = 1;
    matrices.mapped_records = 1;

    iso2gene::Logger logger;
    logger.warn("lf warning");
    iso2gene::write_outputs(config, matrices, logger);

    require_lf_only_file(
        (std::filesystem::path(outdir) / "gene_counts.tsv").string(),
        "gene_counts.tsv"
    );
    require_lf_only_file(
        (std::filesystem::path(outdir) / "gene_tpm.tsv").string(),
        "gene_tpm.tsv"
    );
    require_lf_only_file(
        (std::filesystem::path(outdir) / "gene_length.tsv").string(),
        "gene_length.tsv"
    );
    require_lf_only_file(
        (std::filesystem::path(outdir) / "summary.tsv").string(),
        "summary.tsv"
    );
    require_lf_only_file(
        (std::filesystem::path(outdir) / "warnings.log").string(),
        "warnings.log"
    );

    iso2gene::Logger map_logger;
    const std::string tx2gene_out = temp_output_path("iso2gene_lf_tx2gene.tsv");
    const iso2gene::GtfMapOptions options{
        "tests/data/annotation_basic.gtf",
        tx2gene_out,
        "transcript_id",
        "gene_id"
    };
    (void)iso2gene::make_tx2gene_from_gtf(options, map_logger);
    require_lf_only_file(tx2gene_out, "make-map tx2gene.tsv");
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
        test_version_cli_parse();
        test_make_map_cli_parse();
        test_make_map_basic_gtf();
        test_make_map_gencode_attribute_edges();
        test_make_map_rejects_conflicts_and_bad_attributes();
        test_text_reader_line_endings_and_bom();
        test_external_gzip_fixture();
        test_gzip_optional_header_fields();
        test_gzip_inputs_match_plain();
        test_gzip_concatenated_members_and_errors();
        test_file_outputs_use_lf();
    } catch (const std::exception& error) {
        std::cerr << "TEST FAILED: " << error.what() << '\n';
        return 1;
    }

    std::cout << "all tests passed\n";
    return 0;
}
