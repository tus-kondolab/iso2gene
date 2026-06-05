args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || length(args) > 2) {
  stop("usage: Rscript scripts/benchmark_tximport_vs_iso2gene.R <repo-root> [iterations]", call. = FALSE)
}

suppressPackageStartupMessages(library(tximport))

repo <- normalizePath(args[[1]], mustWork = TRUE)
iterations <- if (length(args) == 2) as.integer(args[[2]]) else 3L
if (is.na(iterations) || iterations < 1) {
  stop("iterations must be a positive integer", call. = FALSE)
}

source_root <- file.path(repo, "build", "tximport_compare")
source_plain <- file.path(source_root, "plain")
source_tx2gene <- file.path(source_root, "tx2gene.tsv")
iso2gene_exe <- file.path(repo, "build-wsl", "iso2gene")

if (!dir.exists(source_plain) || !file.exists(source_tx2gene)) {
  stop("missing build/tximport_compare inputs; run make_tximport_oracle.R first", call. = FALSE)
}
if (!file.exists(iso2gene_exe)) {
  stop("missing build-wsl/iso2gene; build the WSL Release binary first", call. = FALSE)
}

runs <- c("ERR188297", "ERR188088")
formats <- c("kallisto", "salmon", "rsem")
modes <- c(no = "simple-sum", scaledTPM = "scaled-tpm", lengthScaledTPM = "length-scaled-tpm")

bench_root <- file.path(tempdir(), "iso2gene_tximport_bench")
input_root <- file.path(bench_root, "input")
output_root <- file.path(bench_root, "output")
unlink(bench_root, recursive = TRUE, force = TRUE)
dir.create(input_root, recursive = TRUE, showWarnings = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

copy_dir <- function(from, to) {
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(from, to = dirname(to), recursive = TRUE, copy.date = TRUE)
  if (!ok) {
    stop("failed to copy ", from, " to ", to, call. = FALSE)
  }
}

copy_dir(source_plain, file.path(input_root, "plain"))
file.copy(source_tx2gene, file.path(input_root, "tx2gene.tsv"), overwrite = TRUE)

tx2gene_path <- file.path(input_root, "tx2gene.tsv")

quant_path <- function(format, run) {
  if (format == "kallisto") {
    return(file.path(input_root, "plain", "kallisto", run, "abundance.tsv"))
  }
  if (format == "salmon") {
    return(file.path(input_root, "plain", "salmon", run, "quant.sf"))
  }
  if (format == "rsem") {
    return(file.path(input_root, "plain", "rsem", run, paste0(run, ".isoforms.results")))
  }
  stop("unknown format: ", format, call. = FALSE)
}

write_matrix <- function(matrix, path, row_name = "gene_id") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data <- data.frame(row.names(matrix), matrix, check.names = FALSE)
  colnames(data)[1] <- row_name
  write.table(data, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

write_iso_sample_sheet <- function(format) {
  path <- file.path(input_root, paste0("samples_", format, ".tsv"))
  data <- data.frame(
    sample = runs,
    path = vapply(runs, function(run) quant_path(format, run), character(1)),
    stringsAsFactors = FALSE
  )
  write.table(data, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  path
}

sample_sheets <- setNames(vapply(formats, write_iso_sample_sheet, character(1)), formats)

time_expr <- function(expr) {
  gc(FALSE)
  start <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - start
}

run_tximport <- function(format, mode_name, rep_id) {
  files <- setNames(vapply(runs, function(run) quant_path(format, run), character(1)), runs)
  outdir <- file.path(output_root, "tximport", format, mode_name, paste0("rep", rep_id))
  unlink(outdir, recursive = TRUE, force = TRUE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  time_expr({
    tx2gene <- read.delim(tx2gene_path, stringsAsFactors = FALSE)
    txi_args <- list(
      files = files,
      type = format,
      tx2gene = tx2gene,
      countsFromAbundance = mode_name,
      ignoreAfterBar = TRUE,
      dropInfReps = TRUE,
      importer = read.delim
    )
    if (format == "rsem") {
      txi_args$txIn <- TRUE
      txi_args$txOut <- FALSE
    }
    txi <- suppressMessages(do.call(tximport::tximport, txi_args))
    write_matrix(txi$counts, file.path(outdir, "gene_counts.tsv"))
    write_matrix(txi$abundance, file.path(outdir, "gene_tpm.tsv"))
    write_matrix(txi$length, file.path(outdir, "gene_length.tsv"))
  })
}

run_iso2gene <- function(format, mode_name, rep_id) {
  outdir <- file.path(output_root, "iso2gene", format, mode_name, paste0("rep", rep_id))
  unlink(outdir, recursive = TRUE, force = TRUE)
  args <- c(
    "counts",
    "--type", format,
    "--map", tx2gene_path,
    "--sample-sheet", sample_sheets[[format]],
    "--mode", modes[[mode_name]],
    "--outdir", outdir,
    "--ignore-after-bar",
    "--precision", "17"
  )

  time_expr({
    status <- system2(iso2gene_exe, args = args, stdout = FALSE, stderr = FALSE)
    if (!identical(status, 0L)) {
      stop("iso2gene failed with exit status ", status, call. = FALSE)
    }
  })
}

rows <- list()
for (format in formats) {
  for (mode_name in names(modes)) {
    for (rep_id in seq_len(iterations)) {
      rows[[length(rows) + 1]] <- data.frame(
        tool = "tximport",
        format = format,
        mode = mode_name,
        rep = rep_id,
        seconds = run_tximport(format, mode_name, rep_id),
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1]] <- data.frame(
        tool = "iso2gene",
        format = format,
        mode = mode_name,
        rep = rep_id,
        seconds = run_iso2gene(format, mode_name, rep_id),
        stringsAsFactors = FALSE
      )
    }
  }
}

raw <- do.call(rbind, rows)
summary <- aggregate(seconds ~ tool + format + mode, raw, function(x) {
  c(median = median(x), min = min(x), max = max(x))
})
summary <- do.call(data.frame, summary)
colnames(summary) <- c("tool", "format", "mode", "median_seconds", "min_seconds", "max_seconds")

wide <- reshape(
  summary[, c("tool", "format", "mode", "median_seconds")],
  idvar = c("format", "mode"),
  timevar = "tool",
  direction = "wide"
)
colnames(wide) <- sub("^median_seconds\\.", "median_", colnames(wide))
if (all(c("median_tximport", "median_iso2gene") %in% colnames(wide))) {
  wide$speedup_tximport_over_iso2gene <- wide$median_tximport / wide$median_iso2gene
}

dir.create(file.path(repo, "build", "benchmark"), recursive = TRUE, showWarnings = FALSE)
write.table(
  raw,
  file = file.path(repo, "build", "benchmark", "tximport_vs_iso2gene_raw.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  summary,
  file = file.path(repo, "build", "benchmark", "tximport_vs_iso2gene_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  wide,
  file = file.path(repo, "build", "benchmark", "tximport_vs_iso2gene_wide.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

print(summary, row.names = FALSE)
cat("\nMedian speedup ratio, tximport / iso2gene:\n")
print(wide, row.names = FALSE)
