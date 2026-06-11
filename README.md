# iso2gene: compact transcript-to-gene summarization tool

Current version: `1.0.0`

`iso2gene` is a clean-room C/C++ implementation of the core functionality of the **`tximport`** package. iso2gene is provided as a compact command line tool and offers native support for Windows, Linux, and macOS.

In addition to tximport-style summarization, iso2gene can generate
`tx2gene.tsv` directly from GTF annotations. This `make-map`
subcommand is an iso2gene convenience feature and is not part of tximport
itself.

It is designed for workflows that need **`tximport`**-compatible gene-level matrices from a
small standalone command line tool:

```text
kallisto / Salmon / RSEM transcript quantification
+ transcript-to-gene map
-> gene-level counts, TPM, and effective length matrices
-> DESeq2 / PyDESeq2 / edgeR / downstream RNA-seq analysis
```

## Features

- Produces gene-level counts, TPM, and effective length matrices from kallisto
  `abundance.tsv`, Salmon `quant.sf`, and RSEM `isoforms.results`
- Matches tximport `countsFromAbundance="no"`, `"scaledTPM"`, and
  `"lengthScaledTPM"` for supported inputs
- Produces DTU-oriented transcript-level matrices matching tximport
  `countsFromAbundance="dtuScaledTPM"` with `txOut=TRUE`
- Reads plain text and `.gz` quantification files directly
- Generates `tx2gene.tsv` from plain text or gzip-compressed GTF annotations
- Runs as a native command line executable on Windows, Linux, and macOS

## Current Limitations

The following tximport features are not currently implemented in iso2gene:

- kallisto HDF5 input, `abundance.h5`
- tximport input types `sailfish`, `stringtie`, `alevin`, `piscem`, and
  `oarfish`
- tximport custom-column input mode, `type="none"`
- Inferential replicates, such as Salmon Gibbs/bootstrap or kallisto bootstrap
  data

## Install

### Windows

Download the Windows release ZIP from the GitHub Releases page:

https://github.com/tus-kondolab/iso2gene/releases

For version `1.0.0`, download:

```text
iso2gene-1.0.0-windows-x64.zip
```

Extract the ZIP file and run `iso2gene.exe` from PowerShell:

```powershell
.\iso2gene.exe --version
```

### Build From Source

Use these commands to build iso2gene locally instead of using a prebuilt
release binary.

Linux:

```bash
sudo apt update
sudo apt install -y git g++ cmake ninja-build
git clone https://github.com/tus-kondolab/iso2gene.git
cd iso2gene
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

macOS:

```bash
xcode-select --install
brew install git cmake ninja
git clone https://github.com/tus-kondolab/iso2gene.git
cd iso2gene
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

Windows:

Use Developer PowerShell for VS. Visual Studio with the "Desktop development with C++" workload must be installed.

```powershell
git clone https://github.com/tus-kondolab/iso2gene.git
cd iso2gene
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

The commands above build Release binaries.

## Quick Start

Generate a transcript-to-gene map from a GTF annotation:

```bash
iso2gene make-map \
  --gtf gencode.annotation.gtf \
  --out tx2gene.tsv
```

If you already have a two-column `tx2gene.tsv` file, you can use it directly
and skip this step.

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
.\iso2gene.exe counts `
  --type kallisto `
  --map tx2gene.tsv `
  --sample-sheet samples.tsv `
  --mode length-scaled-tpm `
  --outdir out
```

## Supported Input Formats

| `--type` | File | Required columns | Notes |
|---|---|---|---|
| `kallisto` | `abundance.tsv`, `abundance.tsv.gz` | `target_id`, `length`, `eff_length`, `est_counts`, `tpm` | kallisto TSV output |
| `salmon` | `quant.sf`, `quant.sf.gz` | `Name`, `Length`, `EffectiveLength`, `TPM`, `NumReads` | Salmon transcript quantification |
| `rsem` | `isoforms.results`, `isoforms.results.gz` | `transcript_id`, `length`, `effective_length`, `expected_count`, `TPM` | transcript-level RSEM input |

Column names are matched exactly, including case. For example, Salmon uses
`TPM` and `NumReads`, while kallisto uses lowercase `tpm`.

For RSEM input, `effective_length` values smaller than `1` are clamped to `1`
to match tximport's RSEM behavior.

## GZIP Input

iso2gene reads plain text and `.gz` inputs directly. This applies to
quantification files, `tx2gene` maps, sample sheets, and GTF files used by `make-map`.

## Command Line

```text
iso2gene counts --type TYPE --map tx2gene.tsv --sample-sheet samples.tsv --outdir out [options]
iso2gene counts --type TYPE --map tx2gene.tsv --outdir out [options] sample=quant-file ...
iso2gene txout --type TYPE --map tx2gene.tsv --sample-sheet samples.tsv --outdir out [options]
iso2gene txout --type TYPE --map tx2gene.tsv --outdir out [options] sample=quant-file ...
iso2gene make-map --gtf annotation.gtf[.gz] --out tx2gene.tsv [options]
iso2gene --version
```

Required for `counts` and `txout`:

| Option | Description |
|---|---|
| `--map PATH` | Two-column transcript-to-gene TSV |
| `--sample-sheet PATH` | TSV with `sample` and `path` columns, unless direct sample inputs are used |

Options:

| Option | Description |
|---|---|
| `--type TYPE` | `kallisto`, `salmon`, or `rsem` |
| `--mode MODE` | For `counts`: `simple-sum`, `scaled-tpm`, or `length-scaled-tpm`; for `txout`: `dtu-scaled-tpm` |
| `--outdir DIR` | Output directory; default is `out` |
| `--precision N` | Numeric output precision from `1` to `17`; default is `10` |
| `--ignore-version` | Strip transcript suffix after the first dot |
| `--ignore-after-bar` | Strip transcript suffix after the first bar |
| `--help` | Show help |

Required for `make-map`:

| Option | Description |
|---|---|
| `--gtf PATH` | Plain text or gzip-compressed GTF annotation |
| `--out PATH` | Output `tx2gene.tsv` path |

Options for `make-map`:

| Option | Description |
|---|---|
| `--transcript-id-attr NAME` | GTF attribute name for transcript IDs; default is `transcript_id` |
| `--gene-id-attr NAME` | GTF attribute name for gene IDs; default is `gene_id` |
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

### `dtu-scaled-tpm`

This mode is transcript-level and is available through `iso2gene txout`, not
through `iso2gene counts`.

It corresponds to tximport `countsFromAbundance="dtuScaledTPM"` with
`txOut=TRUE`. It is intended for downstream DTU analysis, not for gene-level
DE analysis.

Example:

```bash
iso2gene txout \
  --type salmon \
  --map tx2gene.tsv \
  --sample-sheet samples.tsv \
  --mode dtu-scaled-tpm \
  --outdir out_dtu
```

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

For `iso2gene txout --outdir out_dtu`, iso2gene writes:

```text
out_dtu/transcript_counts.tsv
out_dtu/transcript_tpm.tsv
out_dtu/transcript_length.tsv
out_dtu/transcript_gene.tsv
out_dtu/summary.tsv
out_dtu/warnings.log
```

These matrices use transcripts as rows and samples as columns.

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

`iso2gene make-map` writes transcript IDs as they appear in the GTF. If version
or bar suffix normalization is needed, apply `--ignore-version` or
`--ignore-after-bar` when running `counts`, so the same rule is applied to both
the quantification file and `tx2gene.tsv`.

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

To create this file from a plain text or gzip-compressed GTF annotation:

```bash
iso2gene make-map \
  --gtf annotation.gtf.gz \
  --out tx2gene.tsv
```

`make-map` extracts rows that contain both `transcript_id` and `gene_id`
attributes, collapses duplicate identical transcript-to-gene mappings, and
rejects transcripts that map to multiple genes. GFF3 input is not supported.

## Validation Against tximport

The repository includes scripts for comparing iso2gene with Bioconductor
tximport:

```bash
Rscript scripts/make_tximport_oracle.R /path/to/iso2gene
Rscript scripts/compare_tximport_oracle.R /path/to/iso2gene
```

For transcript-level `dtuScaledTPM` / `txout` validation:

```bash
Rscript scripts/make_txout_oracle.R /path/to/iso2gene
Rscript scripts/compare_txout_oracle.R /path/to/iso2gene oracle_tximport
```

If pytximport is installed, the same iso2gene `txout` output can also be
compared with pytximport:

```bash
python scripts/make_pytximport_txout_oracle.py /path/to/iso2gene
Rscript scripts/compare_txout_oracle.R /path/to/iso2gene oracle_pytximport
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

The oracle script reads the tximportData `.gz` inputs directly and also writes
`build/tximport_compare/samples_gz_<format>.tsv` sample sheets for validating
iso2gene's direct gzip input path.

The validation covers:

```text
kallisto / salmon / rsem
countsFromAbundance = no / scaledTPM / lengthScaledTPM
counts / TPM / effective length matrices
txOut=TRUE with countsFromAbundance = dtuScaledTPM
transcript counts / TPM / effective length matrices
```

With tximportData-derived two-sample fixtures, iso2gene has been checked
against tximport to within `1e-6` absolute and relative tolerance when using
`--precision 17` for output.

The current validation target is Bioconductor tximport `1.38.2`. This matters
for RSEM input because tximport applies an RSEM-specific transcript length
clamp before gene-level summarization.

## Third-Party Code

iso2gene vendors miniz for gzip input support. miniz is MIT licensed; see
`THIRD_PARTY_NOTICES.txt` and `LICENSES/miniz-LICENSE.txt`.

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
