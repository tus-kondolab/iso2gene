args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(
    "usage: Rscript scripts/compare_txout_oracle.R <repo-root> <oracle-name>",
    call. = FALSE
  )
}

repo <- normalizePath(args[[1]], mustWork = TRUE)
oracle_name <- args[[2]]
root <- file.path(repo, "build", "txout_compare")
formats <- c("kallisto", "salmon", "rsem")
matrices <- c(
  counts = "transcript_counts.tsv",
  tpm = "transcript_tpm.tsv",
  length = "transcript_length.tsv"
)

read_matrix <- function(path) {
  data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(data) < 2) {
    stop("matrix has fewer than two columns: ", path, call. = FALSE)
  }
  rownames(data) <- data[[1]]
  data[[1]] <- NULL
  as.matrix(data)
}

compare_one <- function(oracle, observed) {
  if (!setequal(rownames(oracle), rownames(observed))) {
    missing_observed <- setdiff(rownames(oracle), rownames(observed))
    missing_oracle <- setdiff(rownames(observed), rownames(oracle))
    stop(
      "row names differ; missing observed=", length(missing_observed),
      ", missing oracle=", length(missing_oracle),
      call. = FALSE
    )
  }
  if (!setequal(colnames(oracle), colnames(observed))) {
    stop("column names differ", call. = FALSE)
  }

  observed <- observed[rownames(oracle), colnames(oracle), drop = FALSE]
  finite <- is.finite(oracle) & is.finite(observed)
  if (any(is.na(oracle) != is.na(observed))) {
    stop("NA pattern differs", call. = FALSE)
  }
  if (!any(finite)) {
    return(c(max_abs = 0, max_rel = 0))
  }
  diff <- abs(oracle[finite] - observed[finite])
  denom <- pmax(1, abs(oracle[finite]))
  c(max_abs = max(diff), max_rel = max(diff / denom))
}

rows <- list()
failed <- FALSE
for (format in formats) {
  for (matrix_name in names(matrices)) {
    file_name <- matrices[[matrix_name]]
    oracle_path <- file.path(root, oracle_name, format, file_name)
    observed_path <- file.path(root, "iso", format, file_name)

    status <- "OK"
    max_abs <- NA_real_
    max_rel <- NA_real_
    error <- ""
    tryCatch({
      stats <- compare_one(read_matrix(oracle_path), read_matrix(observed_path))
      max_abs <- stats[["max_abs"]]
      max_rel <- stats[["max_rel"]]
      if (max_abs > 1e-6 && max_rel > 1e-6) {
        status <- "FAIL"
        failed <<- TRUE
      }
    }, error = function(e) {
      status <<- "FAIL"
      error <<- conditionMessage(e)
      failed <<- TRUE
    })

    rows[[length(rows) + 1]] <- data.frame(
      oracle = oracle_name,
      format = format,
      matrix = matrix_name,
      status = status,
      max_abs = max_abs,
      max_rel = max_rel,
      error = error,
      stringsAsFactors = FALSE
    )
  }
}

report <- do.call(rbind, rows)
write.table(
  report,
  file = file.path(root, paste0("comparison_", oracle_name, ".tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
print(report, row.names = FALSE)

if (failed) {
  quit(status = 1)
}
