#!/usr/bin/env python3
"""Create pytximport dtu_scaled_tpm transcript-level oracle matrices.

This script is for development validation only. It expects the tximportData-
derived files to be present under tests/data/extdata and the plain extracted
inputs to exist under build/tximport_compare/plain. Run
scripts/make_tximport_oracle.R first if those plain inputs are missing.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pytximport
from pytximport import tximport


RUNS = ["ERR188297", "ERR188088"]
FORMATS = ["kallisto", "salmon", "rsem"]


def plain_paths(repo: Path, data_type: str) -> list[str]:
    root = repo / "build" / "tximport_compare" / "plain" / data_type
    if data_type == "kallisto":
        return [str(root / run / "abundance.tsv") for run in RUNS]
    if data_type == "salmon":
        return [str(root / run / "quant.sf") for run in RUNS]
    if data_type == "rsem":
        return [str(root / run / f"{run}.isoforms.results") for run in RUNS]
    raise ValueError(f"unknown format: {data_type}")


def write_matrix(dataset, variable: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    values = dataset[variable].values
    transcript_ids = [str(value) for value in dataset.coords["transcript_id"].values]
    frame = pd.DataFrame(values, index=transcript_ids, columns=RUNS)
    frame.insert(0, "transcript_id", frame.index)
    frame.to_csv(path, sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("repo_root")
    args = parser.parse_args()

    repo = Path(args.repo_root).resolve()
    out_root = repo / "build" / "txout_compare"
    oracle_root = out_root / "oracle_pytximport"
    tx2gene_path = out_root / "tx2gene.tsv"

    if not tx2gene_path.exists():
        raise SystemExit(
            f"missing {tx2gene_path}; run scripts/make_txout_oracle.R first"
        )

    print(f"pytximport {getattr(pytximport, '__version__', 'unknown')}")
    for data_type in FORMATS:
        files = plain_paths(repo, data_type)
        missing = [path for path in files if not Path(path).exists()]
        if missing:
            raise SystemExit(
                "missing plain input files; run scripts/make_tximport_oracle.R first: "
                + ", ".join(missing)
            )

        dataset = tximport(
            files,
            data_type=data_type,
            transcript_gene_map=str(tx2gene_path),
            counts_from_abundance="dtu_scaled_tpm",
            return_transcript_data=True,
            ignore_after_bar=True,
            ignore_transcript_version=False,
            output_type="xarray",
        )

        format_dir = oracle_root / data_type
        write_matrix(dataset, "counts", format_dir / "transcript_counts.tsv")
        write_matrix(dataset, "abundance", format_dir / "transcript_tpm.tsv")
        write_matrix(dataset, "length", format_dir / "transcript_length.tsv")

    print(f"wrote pytximport txout oracle to {oracle_root}")


if __name__ == "__main__":
    main()
