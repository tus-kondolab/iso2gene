args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: Rscript scripts/make_tximport_oracle.R <repo-root>", call. = FALSE)
}

suppressPackageStartupMessages(library(tximport))

# Validation dataset:
# - Bioconductor tximportData package, extdata directory
# - tximportData contains quantifier outputs for six GEUVADIS Project RNA-seq
#   samples; this script uses the two runs below.
# - The source files are not bundled with iso2gene. Place the tximportData
#   extdata-derived files under tests/data/extdata before running this script.
expected_tximport_version <- "1.38.2"
actual_tximport_version <- as.character(packageVersion("tximport"))
if (actual_tximport_version != expected_tximport_version) {
  warning(
    "This oracle was developed against tximport ",
    expected_tximport_version,
    "; current tximport is ",
    actual_tximport_version,
    ". Review differences before treating failures as iso2gene regressions.",
    call. = FALSE
  )
}

repo <- normalizePath(args[[1]], mustWork = TRUE)
extdata <- file.path(repo, "tests", "data", "extdata")
out_root <- file.path(repo, "build", "tximport_compare")
oracle_root <- file.path(out_root, "oracle")
plain_root <- file.path(out_root, "plain")

runs <- c("ERR188297", "ERR188088")
formats <- c("kallisto", "salmon", "rsem")
modes <- c(no = "no", scaledTPM = "scaledTPM", lengthScaledTPM = "lengthScaledTPM")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(oracle_root, recursive = TRUE, showWarnings = FALSE)
dir.create(plain_root, recursive = TRUE, showWarnings = FALSE)

tx2gene <- read.csv(file.path(extdata, "tx2gene.gencode.v27.csv"), stringsAsFactors = FALSE)
if (ncol(tx2gene) < 2) {
  stop("tx2gene.gencode.v27.csv must have at least two columns", call. = FALSE)
}
tx2gene <- tx2gene[, 1:2]
colnames(tx2gene) <- c("transcript_id", "gene_id")
write.table(
  tx2gene,
  file = file.path(out_root, "tx2gene.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

gz_path <- function(format, run) {
  if (format == "kallisto") {
    return(file.path(extdata, "kallisto", run, "abundance.tsv.gz"))
  }
  if (format == "salmon") {
    return(file.path(extdata, "salmon", run, "quant.sf.gz"))
  }
  if (format == "rsem") {
    return(file.path(extdata, "rsem", run, paste0(run, ".isoforms.results.gz")))
  }
  stop("unknown format: ", format, call. = FALSE)
}

plain_name <- function(format, run) {
  if (format == "kallisto") {
    return("abundance.tsv")
  }
  if (format == "salmon") {
    return("quant.sf")
  }
  if (format == "rsem") {
    return(paste0(run, ".isoforms.results"))
  }
  stop("unknown format: ", format, call. = FALSE)
}

copy_gz_to_plain <- function(src, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  input <- gzfile(src, open = "rt")
  on.exit(close(input), add = TRUE)
  output <- file(dest, open = "wt")
  on.exit(close(output), add = TRUE)
  repeat {
    lines <- readLines(input, n = 10000, warn = FALSE)
    if (!length(lines)) {
      break
    }
    writeLines(lines, output, useBytes = TRUE)
  }
}

write_matrix <- function(matrix, path, row_name = "gene_id") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data <- data.frame(row.names(matrix), matrix, check.names = FALSE)
  colnames(data)[1] <- row_name
  write.table(data, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

for (format in formats) {
  sample_paths <- character(length(runs))
  names(sample_paths) <- runs
  gz_sample_paths <- character(length(runs))
  names(gz_sample_paths) <- runs

  for (run in runs) {
    src <- gz_path(format, run)
    if (!file.exists(src)) {
      stop("missing input file: ", src, call. = FALSE)
    }

    dest <- file.path(plain_root, format, run, plain_name(format, run))
    copy_gz_to_plain(src, dest)
    sample_paths[[run]] <- dest
    gz_sample_paths[[run]] <- src
  }

  sample_sheet <- data.frame(
    sample = runs,
    path = file.path("build", "tximport_compare", "plain", format, runs, basename(sample_paths)),
    stringsAsFactors = FALSE
  )
  write.table(
    sample_sheet,
    file = file.path(out_root, paste0("samples_", format, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )

  gz_sample_sheet <- data.frame(
    sample = runs,
    path = gz_sample_paths,
    stringsAsFactors = FALSE
  )
  write.table(
    gz_sample_sheet,
    file = file.path(out_root, paste0("samples_gz_", format, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )

  for (mode_name in names(modes)) {
    mode <- modes[[mode_name]]
    txi_args <- list(
      files = gz_sample_paths,
      type = format,
      tx2gene = tx2gene,
      countsFromAbundance = mode,
      ignoreAfterBar = TRUE,
      dropInfReps = TRUE
    )
    if (format == "rsem") {
      txi_args$txIn <- TRUE
      txi_args$txOut <- FALSE
    }

    txi <- do.call(tximport::tximport, txi_args)
    mode_dir <- file.path(oracle_root, format, mode_name)
    write_matrix(txi$counts, file.path(mode_dir, "gene_counts.tsv"))
    write_matrix(txi$abundance, file.path(mode_dir, "gene_tpm.tsv"))
    write_matrix(txi$length, file.path(mode_dir, "gene_length.tsv"))
  }
}

cat("wrote tximport oracle and plain inputs to ", out_root, "\n", sep = "")
