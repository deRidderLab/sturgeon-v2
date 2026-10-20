"""
Convert ``modkit extract`` per-read methylation calls to the bed format used by
``sturgeon-v2 single`` / ``sturgeon-v2 multi``.

Ported from the ``inputtobed`` utility of the original Sturgeon package
(https://github.com/UMCUGenetics/sturgeon). The mapping logic is unchanged and
uses the original defaults: margin 25 bp, unmethylated below 0.3, methylated
above 0.7, 5mC code 'm'. Probe locations always come from the model zip.

- One modkit file   -> one bed file (for ``sturgeon-v2 single``)
- A directory       -> one bed file per modkit file (for ``sturgeon-v2 multi``)

Output bed columns (tab separated, with header):
    chrom  chromStart  chromEnd  methylation_call  probe_id
where methylation_call is 1 (methylated) or 0 (unmethylated).
"""

import logging
import os
import zipfile
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

MARGIN = 25
NEG_THRESHOLD = 0.3
POS_THRESHOLD = 0.7
FIVEMC_CODE = "m"

MODEL_PROBES_FILENAME = "probe.bed"
MODKIT_EXTENSIONS = (".txt", ".tsv")
MODKIT_COLUMNS = ["read_id", "chrom", "ref_position", "mod_qual", "mod_code",
                  "canonical_base", "modified_primary_base"]

PROBE_COLUMNS = ["chr", "start", "end", "ID_REF"]
PROBE_COUNT_COLUMNS = ["methylation_calls", "unmethylation_calls", "total_calls"]
BED_COLUMNS = ["chrom", "chromStart", "chromEnd", "methylation_call", "probe_id"]

_CHROM_ORDER = {
    **{str(i): i for i in range(1, 23)},
    "X": 23, "Y": 24, "M": 25, "MT": 25,
}

_PROBE_COLUMN_ALIASES = {
    "chr": ["chr", "chrom", "chromosome"],
    "start": ["start", "chromStart"],
    "end": ["end", "chromEnd"],
    "ID_REF": ["ID_REF", "name", "probe_id"],
}


# ---------------------------------------------------------------------
# Probes from the model
# ---------------------------------------------------------------------
def read_probes_from_model(model_zip: str) -> pd.DataFrame:
    """Read probe.bed from the model zip; return columns chr, start, end, ID_REF.
    This is the same probes file the prediction uses, so probe IDs always match."""
    if not os.path.isfile(model_zip):
        raise ValueError(f"--model not found, given: {model_zip}")
    if not zipfile.is_zipfile(model_zip):
        raise ValueError(f"--model is not a zip file, given: {model_zip}")

    with zipfile.ZipFile(model_zip, "r") as zf:
        members = [
            m for m in zf.namelist()
            if os.path.basename(m) == MODEL_PROBES_FILENAME
            and not m.startswith("__MACOSX/")
        ]
        if len(members) != 1:
            raise ValueError(
                f"Expected exactly one {MODEL_PROBES_FILENAME} in the model zip "
                f"{model_zip}, found {len(members)}."
            )
        with zf.open(members[0]) as handle:
            probes_df = pd.read_csv(handle, sep="\t", header=0, dtype=str)

    rename = {}
    for target, aliases in _PROBE_COLUMN_ALIASES.items():
        found = [c for c in aliases if c in probes_df.columns]
        if found:
            rename[found[0]] = target
        elif target != "end":
            raise ValueError(
                f"{MODEL_PROBES_FILENAME} in {model_zip} has no {target} column "
                f"(accepted names: {aliases}). Columns found: {list(probes_df.columns)}"
            )
    probes_df = probes_df.rename(columns=rename)

    probes_df["start"] = probes_df["start"].astype(np.int64)
    if "end" in probes_df.columns:
        probes_df["end"] = probes_df["end"].astype(np.int64)
    else:
        probes_df["end"] = probes_df["start"] + 1
    probes_df["chr"] = probes_df["chr"].str.replace(r"^chr", "", regex=True)

    logger.info("Using probes file %s from model %s", members[0], model_zip)
    return probes_df[PROBE_COLUMNS].copy()


# ---------------------------------------------------------------------
# Core mapping
# ---------------------------------------------------------------------
def _map_calls_to_probes_chr(probes_df: pd.DataFrame, calls_df: pd.DataFrame) -> pd.DataFrame:
    """For every probe, binarise all calls within +/- MARGIN bp (score >
    POS_THRESHOLD -> +1, score < NEG_THRESHOLD -> -1, else ignored). The
    majority gives one methylated or unmethylated vote; ties give none."""
    probes_df = probes_df[probes_df["start"] > -1]
    probes_df = probes_df.sort_values(["start"]).reset_index(drop=True)
    calls_df = calls_df.sort_values(["reference_pos"]).reset_index(drop=True)

    starts = probes_df["start"].to_numpy() - MARGIN
    ends = starts + MARGIN * 2 + 1

    queries = calls_df["reference_pos"].to_numpy()
    scores = calls_df["score"].to_numpy(dtype=float)

    s = np.searchsorted(queries, starts, "left")
    n = np.searchsorted(queries, ends, "right")

    bin_scores = np.zeros(scores.shape, dtype=np.int64)
    bin_scores[scores > POS_THRESHOLD] = 1
    bin_scores[scores < NEG_THRESHOLD] = -1
    cumulative = np.concatenate(([0], np.cumsum(bin_scores)))

    # sign of the sum of +/-1 votes == median of the non-zero votes
    final_score = np.sign(cumulative[n] - cumulative[s])

    methylation_calls = (final_score == 1).astype(np.int64)
    unmethylation_calls = (final_score == -1).astype(np.int64)
    return probes_df.assign(
        methylation_calls=methylation_calls,
        unmethylation_calls=unmethylation_calls,
        total_calls=methylation_calls + unmethylation_calls,
    )


def _calls_to_bed(calls_df: pd.DataFrame, probes_df: pd.DataFrame) -> pd.DataFrame:
    """Map calls (chr, reference_pos, score) onto all probes and return the bed
    rows, sorted 1..22, X, Y and by position. Probes without a vote are dropped."""
    calls_df = calls_df.assign(
        chr=calls_df["chr"].astype(str).str.replace(r"^chr", "", regex=True)
    )

    per_chr = [
        _map_calls_to_probes_chr(
            probes_df[probes_df["chr"] == chrom],
            calls_df[calls_df["chr"] == chrom],
        )
        for chrom in np.unique(probes_df["chr"])
    ]
    counts = pd.concat(per_chr, ignore_index=True)
    counts = counts[
        (counts["total_calls"] > 0)
        & (counts["methylation_calls"] != counts["unmethylation_calls"])
    ]

    bed_df = pd.DataFrame({
        "chrom": counts["chr"].astype(str),
        "chromStart": counts["start"],
        "chromEnd": counts["end"],
        "methylation_call": (
            counts["methylation_calls"] > counts["unmethylation_calls"]
        ).astype(int),
        "probe_id": counts["ID_REF"],
    })

    key = bed_df["chrom"].map(_CHROM_ORDER).fillna(99)
    return (
        bed_df.assign(_k=key)
        .sort_values(["_k", "chromStart"], kind="stable")
        .drop(columns="_k")
        .reset_index(drop=True)[BED_COLUMNS]
    )


# ---------------------------------------------------------------------
# Modkit
# ---------------------------------------------------------------------
def _is_modkit_file(path: str) -> bool:
    with open(path, "r") as f:
        columns = f.readline().rstrip("\n").split("\t")
    missing = [c for c in MODKIT_COLUMNS if c not in columns]
    if missing:
        logger.error("Skipping %s, not a modkit extract file (missing: %s)", path, missing)
        return False
    return True


def _read_modkit_calls(input_file: str) -> pd.DataFrame:
    modkit_df = pd.read_csv(
        input_file, sep="\t", header=0, usecols=MODKIT_COLUMNS, dtype={"chrom": str}
    )
    modkit_df = modkit_df[
        (modkit_df["mod_code"] == FIVEMC_CODE)
        & (modkit_df["canonical_base"] == "C")
        & (modkit_df["modified_primary_base"] == "C")
    ]
    modkit_df = modkit_df.rename(columns={
        "chrom": "chr", "ref_position": "reference_pos", "mod_qual": "score",
    })
    return modkit_df[(modkit_df["reference_pos"] != -1) & (modkit_df["chr"] != ".")]


def _convert_file(input_file: str, output_file: str, probes_df: pd.DataFrame) -> str:
    logger.info("Processing modkit file: %s", input_file)
    bed_df = _calls_to_bed(_read_modkit_calls(input_file), probes_df)
    logger.info("Total measured array CpG sites: %s", bed_df.shape[0])
    logger.info("Saving bed file to: %s", output_file)
    bed_df.to_csv(output_file, header=True, index=False, sep="\t")
    return output_file


# ---------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------
def modkit_to_bed(input_path: str, output_path: str, model_zip: str) -> List[str]:
    """Convert modkit extract output to sturgeon bed file(s).

    - ``input_path`` is a file: ``output_path`` is the bed file to write
      (if it is an existing directory, ``<name>.bed`` is written inside it).
    - ``input_path`` is a directory: every modkit .txt/.tsv file in it becomes
      ``<output_path>/<name>.bed``, ready for ``sturgeon-v2 multi``.

    Returns the list of written bed files.
    """
    probes_df = read_probes_from_model(model_zip)

    if os.path.isfile(input_path):
        if not _is_modkit_file(input_path):
            raise ValueError(f"{input_path} is not a modkit extract file.")
        if os.path.isdir(output_path):
            output_file = os.path.join(output_path, Path(input_path).stem + ".bed")
        else:
            output_file = output_path
            parent = os.path.dirname(os.path.abspath(output_file))
            os.makedirs(parent, exist_ok=True)
        return [_convert_file(input_path, output_file, probes_df)]

    if os.path.isdir(input_path):
        input_files = [
            os.path.join(input_path, f)
            for f in sorted(os.listdir(input_path))
            if f.endswith(MODKIT_EXTENSIONS)
        ]
        input_files = [f for f in input_files if _is_modkit_file(f)]
        if not input_files:
            raise ValueError(f"No modkit extract files (.txt/.tsv) found in {input_path}")
        logger.info("Found a total of %s modkit files", len(input_files))

        os.makedirs(output_path, exist_ok=True)
        return [
            _convert_file(f, os.path.join(output_path, Path(f).stem + ".bed"), probes_df)
            for f in input_files
        ]

    raise ValueError(f"--input-path must be a file or directory, given: {input_path}")
