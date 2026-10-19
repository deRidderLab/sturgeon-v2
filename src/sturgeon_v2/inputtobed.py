"""
Convert methylation calls (modkit, modkit pileup, megalodon, guppy bam) to the
bed format used by ``sturgeon-v2 single`` / ``sturgeon-v2 multi``.

Ported from the ``inputtobed`` utility of the original Sturgeon package
(https://github.com/UMCUGenetics/sturgeon, sturgeon/cli/inputtobed.py and
sturgeon/callmapping.py). The mapping logic is unchanged; the code has been
updated for pandas >= 2 / numpy >= 2 and the per-probe loop is vectorised.

Only the final bed file is written; all intermediate results stay in memory.

Output bed columns (tab separated, with header):
    chrom  chromStart  chromEnd  methylation_call  probe_id
where methylation_call is 1 (methylated) or 0 (unmethylated).
"""

import gzip
import logging
import os
import shutil
import tempfile
import warnings
import zipfile
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

SOURCES = ["modkit", "modkit_pileup", "megalodon", "guppy"]
MODEL_PROBES_FILENAME = "probe.bed"
MERGED_BED_FILENAME = "merged_probes_methyl_calls.bed"

PROBE_COLUMNS = ["chr", "start", "end", "ID_REF"]
PROBE_COUNT_COLUMNS = ["methylation_calls", "unmethylation_calls", "total_calls"]
BED_COLUMNS = ["chrom", "chromStart", "chromEnd", "methylation_call", "probe_id"]

_CHROM_ORDER = {
    **{str(i): i for i in range(1, 23)},
    "X": 23, "Y": 24, "M": 25, "MT": 25,
}


# ---------------------------------------------------------------------
# Probe file helpers
# ---------------------------------------------------------------------
def extract_probes_from_model(model_zip: str, dest_dir: str) -> str:
    """Extract the probes file (probe.bed) from a model zip into dest_dir.

    This is the same probes file that ``sturgeon-v2 single``/``multi`` use, so
    the probe IDs in the bed file always match the model.
    """
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
        target = os.path.join(dest_dir, MODEL_PROBES_FILENAME)
        with zf.open(members[0]) as src, open(target, "wb") as dst:
            shutil.copyfileobj(src, dst)

    logger.info("Using probes file %s from model %s", members[0], model_zip)
    return target


def _open_text(path: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


# Accepted column names in a probes file, mapped to the names used internally.
_PROBE_COLUMN_ALIASES = {
    "chr": ["chr", "chrom", "chromosome"],
    "start": ["start", "chromStart"],
    "end": ["end", "chromEnd"],
    "ID_REF": ["ID_REF", "name", "probe_id"],
}


def read_probes_file(probes_file: str) -> pd.DataFrame:
    """Read a probes file and return it with columns chr, start, end, ID_REF.

    Accepts the probe.bed from a model zip, the original Sturgeon probes files
    (space separated) and '#'-header samtools-style bed files. Plain or gzipped.
    Chromosome names are stored without the 'chr' prefix (e.g. '1', 'X').
    Any other columns in the file (score, strand, ...) are dropped.
    """
    with _open_text(probes_file) as handle:
        first_line = handle.readline()

    if first_line.startswith("#"):
        column_header = first_line.lstrip("#").rstrip("\n").split("\t")
        probes_df = pd.read_csv(
            probes_file, header=None, sep="\t", comment="#", names=column_header,
            dtype=str,
        )
    else:
        probes_df = pd.read_csv(probes_file, header=0, sep=r"\s+", dtype=str)

    rename = {}
    for target, aliases in _PROBE_COLUMN_ALIASES.items():
        found = [c for c in aliases if c in probes_df.columns]
        if found:
            rename[found[0]] = target
        elif target != "end":
            raise ValueError(
                f"Probes file {probes_file} has no {target} column "
                f"(accepted names: {aliases}). Columns found: {list(probes_df.columns)}"
            )
    probes_df = probes_df.rename(columns=rename)

    probes_df["start"] = probes_df["start"].astype(np.int64)
    if "end" in probes_df.columns:
        probes_df["end"] = probes_df["end"].astype(np.int64)
    else:
        probes_df["end"] = probes_df["start"] + 1

    probes_df["chr"] = probes_df["chr"].str.replace(r"^chr", "", regex=True)
    return probes_df[PROBE_COLUMNS].copy()


def _validate_header(file_path: str, mandatory_columns: List[str], name: str):
    with open(file_path, "r") as f:
        first_line = f.readline()
    column_names = first_line.rstrip("\n").split("\t")
    for mc in mandatory_columns:
        if mc not in column_names:
            return False, f"{mc} column missing in {name} file"
    return True, None


def validate_megalodon_file(megalodon_file: str) -> Tuple[bool, Optional[str]]:
    return _validate_header(
        megalodon_file,
        ["read_id", "chrm", "strand", "pos", "mod_log_prob", "can_log_prob", "mod_base"],
        "megalodon",
    )


def validate_modkit_file(modkit_file: str) -> Tuple[bool, Optional[str]]:
    return _validate_header(
        modkit_file,
        ["read_id", "chrom", "ref_position", "mod_qual", "mod_code",
         "canonical_base", "modified_primary_base"],
        "modkit",
    )


def _collect_files(input_path: str, extension: str) -> List[str]:
    if os.path.isfile(input_path):
        return [input_path]
    if os.path.isdir(input_path):
        return sorted(
            os.path.join(input_path, f)
            for f in os.listdir(input_path)
            if f.endswith(extension)
        )
    raise ValueError(
        f"--input-path must be a directory or file, given: {input_path}"
    )


# ---------------------------------------------------------------------
# Core mapping
# ---------------------------------------------------------------------
def map_methyl_calls_to_probes_chr(
    probes_df: pd.DataFrame,
    methyl_calls_per_read: pd.DataFrame,
    margin: int,
    neg_threshold: float,
    pos_threshold: float,
) -> pd.DataFrame:
    """Map per-read methylation calls to probe locations on one chromosome.

    For every probe, all calls within +/- ``margin`` bp are binarised
    (score > pos_threshold -> +1, score < neg_threshold -> -1, otherwise
    ignored). The probe gets one methylated or unmethylated vote from the
    majority of those calls (the median of the +/-1 values in the original
    implementation); ties give no vote.
    """
    probes_df = probes_df[probes_df["start"] > -1]
    probes_df = probes_df.sort_values(["start"]).reset_index(drop=True)
    methyl_calls_per_read = methyl_calls_per_read.sort_values(
        ["reference_pos"]
    ).reset_index(drop=True)

    starts = probes_df["start"].to_numpy() - margin
    ends = starts + margin * 2 + 1

    queries = methyl_calls_per_read["reference_pos"].to_numpy()
    scores = methyl_calls_per_read["score"].to_numpy(dtype=float)

    s = np.searchsorted(queries, starts, "left")
    n = np.searchsorted(queries, ends, "right")

    bin_scores = np.zeros(scores.shape, dtype=np.int64)
    bin_scores[scores > pos_threshold] = 1
    bin_scores[scores < neg_threshold] = -1
    cumulative = np.concatenate(([0], np.cumsum(bin_scores)))

    # sign of the sum of +/-1 votes == median of the non-zero votes
    final_score = np.sign(cumulative[n] - cumulative[s])

    methylation_calls = probes_df["methylation_calls"].to_numpy().copy()
    unmethylation_calls = probes_df["unmethylation_calls"].to_numpy().copy()
    methylation_calls += (final_score == 1)
    unmethylation_calls += (final_score == -1)

    return probes_df.assign(
        methylation_calls=methylation_calls,
        unmethylation_calls=unmethylation_calls,
        total_calls=methylation_calls + unmethylation_calls,
    )


def _calls_to_probes(
    calls_df: pd.DataFrame,
    probes_file: str,
    margin: int,
    neg_threshold: float,
    pos_threshold: float,
) -> pd.DataFrame:
    """Map a dataframe with columns chr, reference_pos, score onto all probes,
    chromosome by chromosome. Chromosome names may have a 'chr' prefix or not."""
    probes_df = read_probes_file(probes_file)
    for column in PROBE_COUNT_COLUMNS:
        probes_df[column] = 0

    calls_df = calls_df.assign(
        chr=calls_df["chr"].astype(str).str.replace(r"^chr", "", regex=True)
    )

    calls_per_probe = []
    for chrom in np.unique(probes_df["chr"]):
        calls_per_probe_chr = map_methyl_calls_to_probes_chr(
            probes_df=probes_df[probes_df["chr"] == chrom],
            methyl_calls_per_read=calls_df[calls_df["chr"] == chrom],
            margin=margin,
            neg_threshold=neg_threshold,
            pos_threshold=pos_threshold,
        )
        calls_per_probe.append(calls_per_probe_chr)
        logger.debug(
            "Found a total of %s methylation array sites on chromosome %s",
            calls_per_probe_chr["total_calls"].sum(), chrom,
        )

    return pd.concat(calls_per_probe, ignore_index=True)


def _probe_counts_to_bed(calls_per_probe: pd.DataFrame) -> pd.DataFrame:
    """Turn probe call counts into one binary call per probe.
    Probes without calls, or with equal methylated/unmethylated counts, are
    dropped. Rows are sorted 1..22, X, Y, then by position."""
    df = calls_per_probe[
        (calls_per_probe["total_calls"] > 0)
        & (calls_per_probe["methylation_calls"] != calls_per_probe["unmethylation_calls"])
    ]

    bed_df = pd.DataFrame({
        "chrom": df["chr"].astype(str),
        "chromStart": df["start"],
        "chromEnd": df["end"],
        "methylation_call": (
            df["methylation_calls"] > df["unmethylation_calls"]
        ).astype(int),
        "probe_id": df["ID_REF"],
    })

    key = bed_df["chrom"].map(_CHROM_ORDER).fillna(99)
    return (
        bed_df.assign(_k=key)
        .sort_values(["_k", "chromStart"], kind="stable")
        .drop(columns="_k")
        .reset_index(drop=True)[BED_COLUMNS]
    )


def _write_bed(bed_df: pd.DataFrame, output_file: str) -> None:
    logger.info("Total measured array CpG sites: %s", bed_df.shape[0])
    logger.info("Saving bed file to: %s", output_file)
    bed_df.to_csv(output_file, header=True, index=False, sep="\t")


def _files_to_bed(
    input_files: List[str],
    output_path: str,
    file_to_calls,
    **kwargs,
) -> str:
    """Shared driver: compute probe call counts per input file in memory, sum
    them over all files and write only the final bed file."""
    all_calls = []
    for input_file in input_files:
        logger.info("Getting methylation calls from: %s", input_file)
        calls_per_probe = file_to_calls(input_file=input_file, **kwargs)
        if calls_per_probe is None:
            logger.info("File processing gave an error, skipping: %s", input_file)
            continue
        all_calls.append(calls_per_probe)

    if not all_calls:
        raise ValueError("No valid input files were found.")

    merged = (
        pd.concat(all_calls, ignore_index=True)
        .groupby(PROBE_COLUMNS, as_index=False, sort=False)[PROBE_COUNT_COLUMNS]
        .sum()
    )

    bed_output_file = os.path.join(output_path, MERGED_BED_FILENAME)
    _write_bed(_probe_counts_to_bed(merged), bed_output_file)
    return bed_output_file


# ---------------------------------------------------------------------
# Source specific readers
# ---------------------------------------------------------------------
def modkit_file_to_calls(
    input_file: str,
    probes_file: str,
    margin: int = 25,
    neg_threshold: float = 0.3,
    pos_threshold: float = 0.7,
    fivemc_code: str = "m",
    **_,
) -> pd.DataFrame:
    """``modkit extract`` output (per read calls) -> probe call counts."""
    usecols = ["read_id", "chrom", "ref_position", "mod_qual", "mod_code",
               "canonical_base", "modified_primary_base"]
    modkit_df = pd.read_csv(
        input_file, sep="\t", header=0, usecols=usecols, dtype={"chrom": str}
    )

    modkit_df = modkit_df[
        (modkit_df["mod_code"] == fivemc_code)
        & (modkit_df["canonical_base"] == "C")
        & (modkit_df["modified_primary_base"] == "C")
    ]
    modkit_df = modkit_df.rename(columns={
        "chrom": "chr", "ref_position": "reference_pos", "mod_qual": "score",
    })
    modkit_df = modkit_df[
        (modkit_df["reference_pos"] != -1) & (modkit_df["chr"] != ".")
    ]

    logger.info("Processing modkit file: %s", input_file)
    return _calls_to_probes(modkit_df, probes_file, margin, neg_threshold, pos_threshold)


def mega_file_to_calls(
    input_file: str,
    probes_file: str,
    margin: int = 25,
    neg_threshold: float = 0.3,
    pos_threshold: float = 0.7,
    **_,
) -> pd.DataFrame:
    """Megalodon per-read modified base calls -> probe call counts."""
    mega_calls = pd.read_csv(input_file, header=0, sep="\t", dtype={"chrm": str})
    mega_calls = mega_calls.rename(columns={"chrm": "chr", "pos": "reference_pos"})
    mega_calls["score"] = np.exp(mega_calls["mod_log_prob"])

    logger.info("Processing megalodon file: %s", input_file)
    return _calls_to_probes(mega_calls, probes_file, margin, neg_threshold, pos_threshold)


def bam_file_to_calls(
    input_file: str,
    probes_file: str,
    margin: int = 25,
    neg_threshold: float = 0.3,
    pos_threshold: float = 0.7,
    **_,
) -> pd.DataFrame:
    """Guppy bam with MM/ML tags -> probe call counts (needs modbampy)."""
    try:
        from modbampy import ModBam
    except ImportError as e:
        raise ImportError(
            "The 'guppy' source needs modbampy and pysam. Install them with "
            "`pip install sturgeon-v2[bam]`, or (recommended) use modkit instead."
        ) from e

    probes_df = read_probes_file(probes_file)
    results = {"chr": [], "reference_pos": [], "score": []}

    for chrom in np.unique(probes_df["chr"]):
        chromosome = f"chr{chrom}"
        with ModBam(input_file) as bam:
            for read in bam.reads(chromosome, 0, 250000000):
                for pos_mod in read.mod_sites:
                    _, reference_pos, _, _, _, _, _, score = pos_mod
                    if reference_pos == -1:
                        continue
                    results["chr"].append(chromosome)
                    results["reference_pos"].append(reference_pos)
                    results["score"].append((1 + score) / 256)

    logger.info("Processing bam file: %s", input_file)
    return _calls_to_probes(
        pd.DataFrame(results), probes_file, margin, neg_threshold, pos_threshold
    )


def modkit_pileup_file_to_bed(
    input_file: str,
    output_file: str,
    probes_file: str,
    margin: int = 25,
    neg_threshold: float = 0.3,
    pos_threshold: float = 0.7,
    fivemc_code: str = "m",
) -> pd.DataFrame:
    """``modkit pileup`` bedMethyl file -> bed file (written to output_file)."""
    column_names = [
        "chrom", "chromStart", "chromEnd", "mod_code", "score_bed", "strand",
        "thickStart", "thickEnd", "color", "valid_cov", "percent_modified",
        "n_mod", "n_canonical", "n_othermod", "n_delete", "n_fail", "n_diff",
        "n_nocall",
    ]
    modkit_df = pd.read_csv(input_file, sep=r"\s+", header=None, dtype={0: str})

    if modkit_df.shape[1] != len(column_names):
        raise ValueError(
            "Invalid modkit pileup file, number of columns does not match, "
            f"expected {len(column_names)}, got {modkit_df.shape[1]}."
        )

    modkit_df.columns = column_names
    modkit_df = modkit_df[modkit_df["mod_code"] == fivemc_code]
    modkit_df = modkit_df.rename(columns={
        "chrom": "chr", "chromStart": "reference_pos", "percent_modified": "score",
    })
    modkit_df["score"] = modkit_df["score"] / 100
    modkit_df = modkit_df[
        (modkit_df["reference_pos"] != -1) & (modkit_df["chr"] != ".")
    ]

    logger.info("Processing modkit pileup file: %s", input_file)
    calls_per_probe = _calls_to_probes(
        modkit_df, probes_file, margin, neg_threshold, pos_threshold
    )

    # A single file gives each probe at most one vote, so this gives the same
    # result as the original tool.
    bed_df = _probe_counts_to_bed(calls_per_probe)
    _write_bed(bed_df, output_file)
    return bed_df


# ---------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------
def filetobed(
    input_path: str,
    output_path: str,
    source: str,
    model_zip: Optional[str] = None,
    probes_file: Optional[str] = None,
    margin: int = 25,
    neg_threshold: float = 0.3,
    pos_threshold: float = 0.7,
    fivemc_code: str = "m",
) -> str:
    """Convert methylation calls to a sturgeon bed file.

    The probe locations are taken from the model zip (``probe.bed``), or from
    a custom ``probes_file``. Returns the path of the written bed file. For
    ``modkit_pileup``, ``output_path`` is the bed file itself; for all other
    sources it is a directory and the result is
    ``<output_path>/merged_probes_methyl_calls.bed``. No other files are written.
    """
    logger.info("File to bed program, source: %s", source)

    if (model_zip is None) == (probes_file is None):
        raise ValueError("Give either --model or --probes-file (exactly one).")

    tmpdir = None
    try:
        if model_zip is not None:
            tmpdir = tempfile.mkdtemp(prefix="sturgeon_probes_")
            probes_file = extract_probes_from_model(model_zip, tmpdir)
        elif not os.path.exists(probes_file):
            raise ValueError(f"--probes-file not found, given: {probes_file}")

        return _filetobed(
            input_path=input_path, output_path=output_path, source=source,
            probes_file=probes_file, margin=margin, neg_threshold=neg_threshold,
            pos_threshold=pos_threshold, fivemc_code=fivemc_code,
        )
    finally:
        if tmpdir is not None:
            shutil.rmtree(tmpdir, ignore_errors=True)


def _filetobed(
    input_path: str,
    output_path: str,
    source: str,
    probes_file: str,
    margin: int,
    neg_threshold: float,
    pos_threshold: float,
    fivemc_code: str,
) -> str:
    if not 0 <= neg_threshold <= 1:
        raise ValueError(f"--neg-threshold must be between 0 and 1, given: {neg_threshold}")
    if not 0 <= pos_threshold <= 1:
        raise ValueError(f"--pos-threshold must be between 0 and 1, given: {pos_threshold}")
    if pos_threshold <= neg_threshold:
        raise ValueError(
            "--pos-threshold cannot be smaller or equal to --neg-threshold, "
            f"given: {pos_threshold} and {neg_threshold}"
        )
    if margin < 0:
        raise ValueError(f"--margin must be zero or a positive integer, given: {margin}")
    if source not in SOURCES:
        raise ValueError(f"--source must be one of {SOURCES}, given: {source}")

    common = dict(
        probes_file=probes_file,
        margin=margin,
        neg_threshold=neg_threshold,
        pos_threshold=pos_threshold,
    )

    if source == "modkit_pileup":
        if not os.path.isfile(input_path):
            raise ValueError(f"--input-path must be a file for modkit_pileup, given: {input_path}")
        parent = os.path.dirname(os.path.abspath(output_path))
        os.makedirs(parent, exist_ok=True)
        modkit_pileup_file_to_bed(
            input_file=input_path, output_file=output_path,
            fivemc_code=fivemc_code, **common,
        )
        return output_path

    os.makedirs(output_path, exist_ok=True)

    if source == "modkit":
        candidates = _collect_files(input_path, ".txt")
        input_files = []
        for f in candidates:
            ok, msg = validate_modkit_file(f)
            if ok:
                input_files.append(f)
            else:
                logger.error("File %s did not pass modkit validation: %s", f, msg)
        logger.info("Found a total of %s modkit files", len(input_files))
        return _files_to_bed(
            input_files, output_path, modkit_file_to_calls,
            fivemc_code=fivemc_code, **common,
        )

    if source == "megalodon":
        warnings.warn(
            "Using megalodon is NOT recommended; megalodon is deprecated by ONT. "
            "Please consider modkit (https://github.com/nanoporetech/modkit)."
        )
        candidates = _collect_files(input_path, ".txt")
        input_files = []
        for f in candidates:
            ok, msg = validate_megalodon_file(f)
            if ok:
                input_files.append(f)
            else:
                logger.error("File %s did not pass megalodon validation: %s", f, msg)
        logger.info("Found a total of %s megalodon files", len(input_files))
        return _files_to_bed(input_files, output_path, mega_file_to_calls, **common)

    # source == "guppy"
    warnings.warn(
        "Using guppy bam files is NOT recommended; modbampy is deprecated and "
        "may give WRONG results with newer guppy versions. Please consider "
        "modkit (https://github.com/nanoporetech/modkit)."
    )
    input_files = _collect_files(input_path, ".bam")
    try:
        import pysam
    except ImportError as e:
        raise ImportError(
            "The 'guppy' source needs pysam and modbampy. Install them with "
            "`pip install sturgeon-v2[bam]`, or (recommended) use modkit instead."
        ) from e
    for bam_file in input_files:
        if not os.path.exists(bam_file + ".bai"):
            logger.info("Generating index file: %s.bai", bam_file)
            pysam.index(bam_file)
    logger.info("Found a total of %s bam files", len(input_files))
    return _files_to_bed(input_files, output_path, bam_file_to_calls, **common)
