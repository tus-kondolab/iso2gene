args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: Rscript scripts/make_txout_oracle.R <repo-root>", call. = FALSE)
}

suppressPackageStartupMessages(library(tximport))

# Validation dataset:
# - Bioconductor tximportData package, extdata directory
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
out_root <- file.path(repo, "build", "txout_compare")
oracle_root <- file.path(out_root, "oracle_tximport")

runs <- c("ERR188297", "ERR188088")
formats <- c("kallisto", "salmon", "rsem")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(oracle_root, recursive = TRUE, showWarnings = FALSE)

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

write_matrix <- function(matrix, path, row_name = "transcript_id") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data <- data.frame(row.names(matrix), matrix, check.names = FALSE)
  colnames(data)[1] <- row_name
  write.table(data, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

for (format in formats) {
  files <- character(length(runs))
  names(files) <- runs
  for (run in runs) {
    files[[run]] <- gz_path(format, run)
    if (!file.exists(files[[run]])) {
      stop("missing input file: ", files[[run]], call. = FALSE)
    }
  }

  txi_args <- list(
    files = files,
    type = format,
    tx2gene = tx2gene,
    txOut = TRUE,
    countsFromAbundance = "dtuScaledTPM",
    ignoreAfterBar = TRUE,
    dropInfReps = TRUE
  )
  if (format == "rsem") {
    txi_args$txIn <- TRUE
  }

  txi <- do.call(tximport::tximport, txi_args)
  format_dir <- file.path(oracle_root, format)
  write_matrix(txi$counts, file.path(format_dir, "transcript_counts.tsv"))
  write_matrix(txi$abundance, file.path(format_dir, "transcript_tpm.tsv"))
  write_matrix(txi$length, file.path(format_dir, "transcript_length.tsv"))

  tx_gene <- tx2gene[match(row.names(txi$counts), tx2gene$transcript_id), , drop = FALSE]
  write.table(
    tx_gene,
    file = file.path(format_dir, "transcript_gene.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

cat("wrote txout tximport oracle to ", out_root, "\n", sep = "")
