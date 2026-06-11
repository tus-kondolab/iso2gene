#include "iso2gene/write_matrix.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "iso2gene/error.hpp"

namespace iso2gene {

namespace {

std::string path_join(const std::string& dir, const std::string& filename) {
    return (std::filesystem::path(dir) / filename).string();
}

void write_summary(const Config& config, const GeneMatrices& matrices, const std::string& path) {
    std::ofstream out(path, std::ios::out | std::ios::binary);
    if (!out) {
        throw Iso2GeneError(ExitCode::io_error, "failed to open summary for writing: " + path);
    }

    out << "metric\tvalue\n";
    out << "mode\t" << mode_name(config.mode) << '\n';
    out << "input_type\t" << config.input_type << '\n';
    out << "samples\t" << matrices.sample_names.size() << '\n';
    out << "genes\t" << matrices.gene_ids.size() << '\n';
    out << "total_transcript_records\t" << matrices.total_records << '\n';
    out << "mapped_transcript_records\t" << matrices.mapped_records << '\n';
    out << "unmapped_transcript_records\t" << matrices.unmapped_records << '\n';
}

void write_transcript_summary(
    const Config& config,
    const TranscriptMatrices& matrices,
    const std::string& path
) {
    std::ofstream out(path, std::ios::out | std::ios::binary);
    if (!out) {
        throw Iso2GeneError(ExitCode::io_error, "failed to open summary for writing: " + path);
    }

    out << "metric\tvalue\n";
    out << "mode\t" << transcript_mode_name(config.transcript_mode) << '\n';
    out << "input_type\t" << config.input_type << '\n';
    out << "samples\t" << matrices.sample_names.size() << '\n';
    out << "transcripts\t" << matrices.transcript_ids.size() << '\n';
    out << "total_transcript_records\t" << matrices.total_records << '\n';
    out << "mapped_transcript_records\t" << matrices.mapped_records << '\n';
    out << "unmapped_transcript_records\t" << matrices.unmapped_records << '\n';
}

void write_transcript_gene_tsv(
    const std::string& path,
    const std::vector<std::string>& transcript_ids,
    const std::vector<std::string>& gene_ids
) {
    if (transcript_ids.size() != gene_ids.size()) {
        throw Iso2GeneError(
            ExitCode::internal_error,
            "transcript and gene label vectors have different sizes"
        );
    }

    std::ofstream out(path, std::ios::out | std::ios::binary);
    if (!out) {
        throw Iso2GeneError(ExitCode::io_error, "failed to open transcript-gene map for writing: " + path);
    }

    out << "transcript_id\tgene_id\n";
    for (std::size_t row = 0; row < transcript_ids.size(); ++row) {
        out << transcript_ids[row] << '\t' << gene_ids[row] << '\n';
    }
}

} // namespace

void write_matrix_tsv(
    const std::string& path,
    const std::string& row_name,
    const std::vector<std::string>& row_ids,
    const std::vector<std::string>& col_names,
    const Matrix<double>& matrix,
    int precision
) {
    if (row_ids.size() != matrix.rows() || col_names.size() != matrix.cols()) {
        throw Iso2GeneError(
            ExitCode::internal_error,
            "matrix dimensions do not match row or column labels"
        );
    }

    std::ofstream out(path, std::ios::out | std::ios::binary);
    if (!out) {
        throw Iso2GeneError(ExitCode::io_error, "failed to open matrix for writing: " + path);
    }

    out << std::setprecision(precision);
    out << row_name;
    for (const std::string& col : col_names) {
        out << '\t' << col;
    }
    out << '\n';

    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        out << row_ids[row];
        for (std::size_t col = 0; col < matrix.cols(); ++col) {
            out << '\t' << matrix(row, col);
        }
        out << '\n';
    }
}

void write_outputs(
    const Config& config,
    const GeneMatrices& matrices,
    const Logger& logger
) {
    std::error_code ec;
    std::filesystem::create_directories(config.outdir, ec);
    if (ec) {
        throw Iso2GeneError(
            ExitCode::io_error,
            "failed to create output directory '" + config.outdir + "': " + ec.message()
        );
    }

    write_matrix_tsv(
        path_join(config.outdir, "gene_counts.tsv"),
        "gene_id",
        matrices.gene_ids,
        matrices.sample_names,
        matrices.counts,
        config.precision
    );
    write_matrix_tsv(
        path_join(config.outdir, "gene_tpm.tsv"),
        "gene_id",
        matrices.gene_ids,
        matrices.sample_names,
        matrices.tpm,
        config.precision
    );
    write_matrix_tsv(
        path_join(config.outdir, "gene_length.tsv"),
        "gene_id",
        matrices.gene_ids,
        matrices.sample_names,
        matrices.length,
        config.precision
    );
    write_summary(config, matrices, path_join(config.outdir, "summary.tsv"));
    logger.write_warnings(path_join(config.outdir, "warnings.log"));
}

void write_transcript_outputs(
    const Config& config,
    const TranscriptMatrices& matrices,
    const Logger& logger
) {
    std::error_code ec;
    std::filesystem::create_directories(config.outdir, ec);
    if (ec) {
        throw Iso2GeneError(
            ExitCode::io_error,
            "failed to create output directory '" + config.outdir + "': " + ec.message()
        );
    }

    write_matrix_tsv(
        path_join(config.outdir, "transcript_counts.tsv"),
        "transcript_id",
        matrices.transcript_ids,
        matrices.sample_names,
        matrices.counts,
        config.precision
    );
    write_matrix_tsv(
        path_join(config.outdir, "transcript_tpm.tsv"),
        "transcript_id",
        matrices.transcript_ids,
        matrices.sample_names,
        matrices.tpm,
        config.precision
    );
    write_matrix_tsv(
        path_join(config.outdir, "transcript_length.tsv"),
        "transcript_id",
        matrices.transcript_ids,
        matrices.sample_names,
        matrices.length,
        config.precision
    );
    write_transcript_gene_tsv(
        path_join(config.outdir, "transcript_gene.tsv"),
        matrices.transcript_ids,
        matrices.gene_ids
    );
    write_transcript_summary(config, matrices, path_join(config.outdir, "summary.tsv"));
    logger.write_warnings(path_join(config.outdir, "warnings.log"));
}

} // namespace iso2gene
