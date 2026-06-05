# iso2gene

`iso2gene` is a clean-room C/C++ implementation of the tximport algorighm. It is a small
dependency-free command line tool that summarizes transcript- or isoform-level
quantification files to gene-level matrices.

It is designed for workflows that need tximport-like output without requiring
R at runtime:

```text
kallisto / Salmon / RSEM transcript quantification
+ transcript-to-gene map
-> gene-level counts, TPM, and effective length matrices
-> DESeq2 / PyDESeq2 / edgeR / downstream RNA-seq analysis
```

The implementation is independent from the tximport source code: iso2gene does
not copy, translate, or port tximport internals. Its calculation behavior is
validated against Bioconductor tximport using tximportData-derived fixtures.

## Features

- Supports kallisto `abundance.tsv`
- Supports Salmon `quant.sf`
- Supports RSEM `isoforms.results`
- Produces gene-level counts, TPM, and effective length matrices
- Implements `simple-sum`, `scaled-tpm`, and `length-scaled-tpm`
- Matches tximport `countsFromAbundance="no"`, `"scaledTPM"`, and
  `"lengthScaledTPM"` for supported inputs
- Runs as a native C++17 executable with no external runtime dependencies
- Builds with MSVC, GCC, Clang, and Apple Clang

## Current Limitations

- Input files must be plain text. `.gz` files are not supported yet.
- `abundance.h5` is not supported yet.
- Salmon/RSEM inferential replicates are not imported.
- RSEM `genes.results` is not supported yet.
- GTF/GFF3 parsing and tx2gene generation are not implemented yet.
- `--map` is required for all input types, including RSEM.

## Install

Requirements:

- Git
- CMake 3.20 or newer
- A C++17 compiler

Clone the repository:

```bash
git clone https://github.com/tus-kondolab/iso2gene.git
cd iso2gene
```

Windows with Visual Studio generator:

```powershell
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

This uses the Visual Studio/MSBuild generator and is Windows-specific.

Windows/Linux/macOS with Ninja:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

The executable is typically written to one of:

```text
build/Release/iso2gene.exe
build/iso2gene
```

depending on the generator and platform. The commands above build Release
binaries.

## Quick Start

Create a transcript-to-gene map:

```text
transcript_id	gene_id
ENST00000456328.2	ENSG00000223972.5
ENST00000450305.2	ENSG00000223972.5
```

Create a sample sheet:

```text
sample	path
control_1	kallisto/control_1/abundance.tsv
control_2	kallisto/control_2/abundance.tsv
treated_1	kallisto/treated_1/abundance.tsv
treated_2	kallisto/treated_2/abundance.tsv
```

You may add optional metadata columns such as `condition`, `batch`, or
`replicate`, but iso2gene ignores them when generating gene-level matrices.

Run iso2gene:

```bash
iso2gene counts \
  --type kallisto \
  --map tx2gene.tsv \
  --sample-sheet samples.tsv \
  --mode length-scaled-tpm \
  --outdir out
```

On Windows PowerShell:

```powershell
.\build\Release\iso2gene.exe counts `
  --type kallisto `
  --map tx2gene.tsv `
  --sample-sheet samples.tsv `
  --mode length-scaled-tpm `
  --outdir out
```

## Supported Input Formats

| `--type` | File | Required columns | Notes |
|---|---|---|---|
| `kallisto` | `abundance.tsv` | `target_id`, `length`, `eff_length`, `est_counts`, `tpm` | kallisto TSV output |
| `salmon` | `quant.sf` | `Name`, `Length`, `EffectiveLength`, `TPM`, `NumReads` | Salmon transcript quantification |
| `rsem` | `isoforms.results` | `transcript_id`, `length`, `effective_length`, `expected_count`, `TPM` | transcript-level RSEM input |

Column names are matched exactly, including case. For example, Salmon uses
`TPM` and `NumReads`, while kallisto uses lowercase `tpm`.

For RSEM input, `effective_length` values smaller than `1` are clamped to `1`
to match tximport's RSEM behavior.

## Command Line

```text
iso2gene counts --type TYPE --map tx2gene.tsv --sample-sheet samples.tsv --outdir out [options]
iso2gene counts --type TYPE --map tx2gene.tsv --outdir out [options] sample=quant-file ...
```

Required for `counts`:

| Option | Description |
|---|---|
| `--map PATH` | Two-column transcript-to-gene TSV |
| `--sample-sheet PATH` | TSV with `sample` and `path` columns, unless direct sample inputs are used |

Options:

| Option | Description |
|---|---|
| `--type TYPE` | `kallisto`, `salmon`, or `rsem` |
| `--mode MODE` | `simple-sum`, `scaled-tpm`, or `length-scaled-tpm` |
| `--outdir DIR` | Output directory; default is `out` |
| `--precision N` | Numeric output precision from `1` to `17`; default is `10` |
| `--ignore-version` | Strip transcript suffix after the first dot |
| `--ignore-after-bar` | Strip transcript suffix after the first bar |
| `--help` | Show help |

## Sample Sheet

The recommended way to pass samples is a TSV file:

```text
sample	path
ERR188297	salmon/ERR188297/quant.sf
ERR188088	salmon/ERR188088/quant.sf
```

The required columns are:

- `sample`: sample name used as the output matrix column name
- `path`: path to that sample's quantification file

Additional columns, such as `condition`, are allowed but ignored by iso2gene.

## Direct Inputs

For quick small runs, samples can be passed directly as `sample=path` arguments:

```bash
iso2gene counts \
  --type salmon \
  --map tx2gene.tsv \
  --mode length-scaled-tpm \
  --outdir out \
  control_1=salmon/control_1/quant.sf \
  control_2=salmon/control_2/quant.sf
```

Here `control_1` and `control_2` are sample names, not shell variables. They
become the output matrix column names.

If a path contains spaces, quote the whole argument:

```bash
"control_1=salmon/control 1/quant.sf"
```

## Count Modes

### `simple-sum`

Sums transcript estimated counts by gene:

```text
gene_count[g, s] = sum transcript_count[t, s]
```

This corresponds to tximport `countsFromAbundance="no"` after gene-level
summarization.

### `scaled-tpm`

Uses gene-level TPM and rescales each sample so that the total count matches
the mapped original gene count total.

This corresponds to tximport `countsFromAbundance="scaledTPM"`.

### `length-scaled-tpm`

Uses gene-level TPM multiplied by the sample-averaged gene effective length,
then rescales each sample to the mapped original gene count total.

This corresponds to tximport `countsFromAbundance="lengthScaledTPM"` and is
the recommended mode for many gene-level differential expression workflows.

## Outputs

For `--outdir out`, iso2gene writes:

```text
out/gene_counts.tsv
out/gene_tpm.tsv
out/gene_length.tsv
out/summary.tsv
out/warnings.log
```

Matrix files use genes as rows and samples as columns:

```text
gene_id	control_1	control_2
ENSG000001	100.25	98.75
ENSG000002	0	3.5
```

Counts are not rounded by default. Use downstream tools' recommended rounding
or integer handling policy when necessary.

## Transcript ID Normalization

Some quantification files contain transcript IDs with suffixes that are not
present in `tx2gene.tsv`.

Use `--ignore-version` to strip the suffix after the first dot:

```text
ENST00000456328.2 -> ENST00000456328
```

Use `--ignore-after-bar` to strip the suffix after the first bar:

```text
ENST00000456328.2|ENSG00000223972.5|... -> ENST00000456328.2
```

These transformations are applied to both the quantification file IDs and the
tx2gene transcript IDs.

## tx2gene File

`--map` expects a two-column TSV:

```text
transcript_id	gene_id
ENST00000456328.2	ENSG00000223972.5
ENST00000450305.2	ENSG00000223972.5
```

A header is allowed. If no recognized header is present, the first two columns
are treated as transcript ID and gene ID.

The same transcript may not map to multiple genes. Duplicate identical rows
are ignored with a warning.

## Validation Against tximport

The repository includes scripts for comparing iso2gene with Bioconductor
tximport:

```bash
Rscript scripts/make_tximport_oracle.R /path/to/iso2gene
Rscript scripts/compare_tximport_oracle.R /path/to/iso2gene
```

The validation data is not bundled with iso2gene. The comparison was performed
with files from the Bioconductor `tximportData` package, specifically the
package `extdata` directory containing output from multiple transcript
abundance quantifiers on six GEUVADIS Project RNA-seq samples.

The iso2gene validation used this two-sample subset:

| run | population | center | sample | experiment |
|---|---|---|---|---|
| `ERR188297` | `TSI` | `UNIGE` | `ERS185497` | `ERX163094` |
| `ERR188088` | `TSI` | `UNIGE` | `ERS185242` | `ERX162972` |

For these two runs, the validation used:

```text
tximportData/extdata/tx2gene.gencode.v27.csv
tximportData/extdata/kallisto/<run>/abundance.tsv.gz
tximportData/extdata/salmon/<run>/quant.sf.gz
tximportData/extdata/rsem/<run>/<run>.isoforms.results.gz
```

To run the validation scripts locally, provide those tximportData-derived files
under `tests/data/extdata` with the same directory layout. The external data
directory is not required to build or run iso2gene.

The validation covers:

```text
kallisto / salmon / rsem
countsFromAbundance = no / scaledTPM / lengthScaledTPM
counts / TPM / effective length matrices
```

With tximportData-derived two-sample fixtures, iso2gene has been checked
against tximport to within `1e-6` absolute and relative tolerance when using
`--precision 17` for output.

The current validation target is Bioconductor tximport `1.38.2`. This matters
for RSEM input because tximport applies an RSEM-specific transcript length
clamp before gene-level summarization.

## Development

Run the C++ tests:

```bash
cmake --build build
ctest --test-dir build -C Debug --output-on-failure
```

For single-config generators such as Ninja, omit `-C Debug`:

```bash
ctest --test-dir build --output-on-failure
```

## License

iso2gene is released under the MIT License. See [LICENSE](LICENSE).
