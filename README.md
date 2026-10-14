# Sturgeon V2

**Intraoperative CNS tumour classification from nanopore methylation data.**

Sturgeon V2 is the successor of [Sturgeon](https://github.com/UMCUGenetics/sturgeon), a neural network classifier for central nervous system (CNS) tumours based on sparse DNA methylation profiles, such as those that can be obtained with intra-operative nanopore sequencing. This repository provides a pip-installable package to run Sturgeon V2 in a clinical setting, from preprocessing of methylation calls through to prediction.

## Table of contents
- [Installation](#installation)
- [Models](#models)
- [Workflow overview](#workflow-overview)
- [Step 1: Preprocessing](#step-1-preprocessing)
- [Step 2: Prediction](#step-2-prediction)
- [Command-line reference](#command-line-reference)
- [Output](#output)

## Installation

### Requirements

- Python >= 3.10,<3.15
- Linux or macOS <!-- TODO: confirm supported platforms -->
- [modkit](https://github.com/nanoporetech/modkit) (recommended), for extracting methylation calls from BAM files

The original Sturgeon package is **not** required: preprocessing (`inputtobed`) is included in Sturgeon V2.

### Install Sturgeon V2

We recommend installing into a fresh virtual environment:

```bash
python3 -m venv venv
source venv/bin/activate
python3 -m pip install --upgrade pip
```

Install directly from GitHub over SSH (requires an SSH key linked to your GitHub account):

```bash
pip install git+ssh://git@github.com/UMCUGenetics/sturgeon-v2.git
```

Alternatively, install over HTTPS:

```bash
pip install git+https://github.com/UMCUGenetics/sturgeon-v2.git
```

Or install from a local clone:

```bash
git clone https://github.com/UMCUGenetics/sturgeon-v2.git
cd sturgeon-v2
pip install .
```

Verify the installation:

```bash
sturgeon-v2 --help
```

## Models

Models are distributed as zip files and are passed to `sturgeon-v2` with the `-m/--model` option. Download a model and store it in a location of your choice; it does not need to be unpacked.

| Model | Status | Download |
|---|---|---|
| `cns-v2.4b` | **Current, recommended** | [zip file](https://drive.google.com/file/d/18cHuUTN_7TyOZzcbaIGeVrugoYCHNkOE/view?usp=share_link) |
| `cns-v2.4a` | Deprecated | [zip file](https://drive.google.com/file/d/1opEUbGWtUCONWhlfrv9NfSuGM8pXj-D8/view?usp=share_link) |

<!-- TODO: describe the classification scheme of cns-v2.4b (number of classes, control classes, merged classes) and what changed between v2.4a and v2.4b -->

> **Note:** `cns-v2.4a` is kept available only for reproducibility of earlier results. Use `cns-v2.4b` for all new analyses.

## Workflow overview

```
Nanopore BAM (with 5mC calls, aligned to T2T-CHM13v2.0)
        │
        │  modkit adjust-mods / modkit extract        (Step 1a)
        ▼
modkit output (.txt)
        │
        │  sturgeon-v2 inputtobed                     (Step 1b)
        ▼
Methylation calls on array probe sites (.bed)
        │
        │  sturgeon-v2 single / multi                 (Step 2)
        ▼
Class scores per sample
```

## Step 1: Preprocessing

Preprocessing is done with `sturgeon-v2 inputtobed`. This command is ported from the `inputtobed` utility of the original [Sturgeon](https://github.com/UMCUGenetics/sturgeon) and produces identical output, so `.bed` files previously created with `sturgeon inputtobed` can still be used with Sturgeon V2.

> ⚠️ **Reference genome: only T2T-CHM13v2.0 is supported.**
> Sturgeon V2 does **not** support `hg38`. Reads must be aligned to the [T2T-CHM13v2.0 reference genome](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz), which is the default for `inputtobed`. Do **not** pass `--reference-genome hg38`. If your data is aligned to `hg38`, re-align it to T2T-CHM13v2.0 before preprocessing. Incorrectly aligned data will give wrong results without raising an error.

### 1a. Extract methylation calls with modkit

Starting from a BAM file with modified base calls (from Dorado or Guppy), aligned to T2T-CHM13v2.0:

**(Optional) Combine 5hmC with 5mC.** If your BAM contains both 5mC and 5hmC calls, merge them into a single 5mC score. Both modifications look alike on bisulfite-based methylation arrays, which the reference data was generated with, so combining them best matches the training data.

```bash
modkit adjust-mods --convert h m input.bam adjusted.bam
```

**Extract per-read calls:**

```bash
mkdir -p sample1/modkit
modkit extract full adjusted.bam sample1/modkit/modkit_output.txt
```

Older versions of modkit use `modkit extract` without `full`. Consult the [modkit documentation](https://nanoporetech.github.io/modkit/quick_start.html) if the command has changed.

Alternatively, `modkit pileup` can be used to produce a single file of aggregated calls per genomic position.

### 1b. Convert to bed with `sturgeon-v2 inputtobed`

For `modkit extract` output, pass the directory containing the `.txt` file(s):

```bash
sturgeon-v2 inputtobed -i sample1/modkit -o sample1/bed -s modkit
```

This writes `sample1/bed/merged_probes_methyl_calls.bed`, which is the input for [Step 2](#step-2-prediction). The output directory also contains intermediate files: one `<name>_probes_methyl_calls.txt` per input file and `merged_probes_methyl_calls.txt` with the summed call counts per probe.

For `modkit pileup` output, pass a single pileup file; `-o` is then the output bed file itself:

```bash
sturgeon-v2 inputtobed -i sample1/pileup.bed -o sample1/sample1.bed -s modkit_pileup
```

An intermediate file with the same name plus `.tmp` is written next to it.

> **One sample per run.** `inputtobed` merges all methylation calls from all `.txt` files in the input directory into a single output. If you have multiple samples, process each sample separately into its own output directory.

> **Re-running in the same output directory.** For each input file, `inputtobed` skips processing if `<name>_probes_methyl_calls.txt` already exists in the output directory and reuses it. If input files have changed (for example, during live sequencing, or after changing `--margin` or the thresholds), use a new output directory or delete the old intermediate files first.

### Deprecated input sources

`inputtobed` also accepts Megalodon output (`-s megalodon`) and Guppy BAM files (`-s guppy`, requires the [optional `bam` extra](#optional-support-for-guppy-bam-input)). Megalodon is deprecated by Oxford Nanopore, and Guppy BAM parsing relies on `modbampy`, which is deprecated and may give wrong results with newer Guppy versions. Use modkit whenever possible.

For Guppy BAM files, do not use the `.bai` index files produced by Guppy, as they are not compatible. `inputtobed` creates compatible index files itself when none are present.

## Step 2: Prediction

Sturgeon V2 has two prediction modes:

- `single` predicts on one `.bed` or `.csv` file.
- `multi` predicts on every `.bed` or `.csv` file in a directory, treating each file as an independent sample.

### Predict a single sample

```bash
sturgeon-v2 single \
    -i sample1/bed/merged_probes_methyl_calls.bed \
    -m models/cns-v2.4b.zip \
    -f bed \
    -o results/sample1_prediction.csv
```

### Predict a directory of samples

Collect one bed file per sample in a directory. The file name (without extension) is used as the sample name in the output, so rename the `merged_probes_methyl_calls.bed` files:

```bash
mkdir -p bed_files
cp sample1/bed/merged_probes_methyl_calls.bed bed_files/sample1.bed
cp sample2/bed/merged_probes_methyl_calls.bed bed_files/sample2.bed

sturgeon-v2 multi \
    -i bed_files/ \
    -m models/cns-v2.4b.zip \
    -f bed \
    -o results/predictions.tsv
```

In `multi` mode, `-o` is a single tab-separated output file covering all samples.

### Simulating lower coverage (subsampling)

The `-c` and `-s` options randomly subsample the available CpG sites. This is useful for benchmarking how the classifier behaves at lower sequencing depth, for example to mimic an early timepoint during surgery. Sampling is seeded, so results are reproducible.

The following example draws 10 random subsets of 1,000 CpG sites and predicts on each:

```bash
sturgeon-v2 single \
    -i sample1/bed/merged_probes_methyl_calls.bed \
    -m models/cns-v2.4b.zip \
    -f bed \
    -c 1000 \
    -s 10 \
    -o results/sample1_1000cpgs.csv
```

If `-c` is larger than the number of available CpG sites, all sites are used. For routine clinical use, leave `-c` unset so that all available CpG sites are used.

## Command-line reference

`sturgeon-v2` has three subcommands:

```
sturgeon-v2 inputtobed -i INPUT_PATH -o OUTPUT_PATH -s SOURCE [options]
sturgeon-v2 single     -i INPUT_FILE -m MODEL -o OUTPUT_FILE -f FILE_TYPE [-c NUM_CPGS] [-s NUM_SAMPLINGS]
sturgeon-v2 multi      -i INPUT_DIR  -m MODEL -o OUTPUT_FILE -f FILE_TYPE [-c NUM_CPGS] [-s NUM_SAMPLINGS]
```

Run `sturgeon-v2 <subcommand> --help` for the full list of options.

### `inputtobed`

| Option | Description | Default |
|---|---|---|
| `-i`, `--input-path` | Input file or directory. For `modkit`/`megalodon`, all `.txt` files in the directory are used; for `guppy`, all `.bam` files. For `modkit_pileup`, a single file. | — (required) |
| `-o`, `--output-path` | Output directory. For `modkit_pileup`, the output `.bed` file. | — (required) |
| `-s`, `--source` | Tool that produced the input: `modkit`, `modkit_pileup`, `megalodon` or `guppy`. | — (required) |
| `--reference-genome` | Reference genome used for alignment. Only `chm13v2` is supported by Sturgeon V2 (see [Step 1](#step-1-preprocessing)). | `chm13v2` |
| `--probes-file` | Custom probe location file, instead of the bundled one. Cannot be combined with `--reference-genome`. | Bundled CHM13v2 probes |
| `--margin` | Distance (bp) around each probe within which methylation calls are considered. | `25` |
| `--neg-threshold` | Calls with a score below this value are counted as unmethylated. | `0.3` |
| `--pos-threshold` | Calls with a score above this value are counted as methylated. | `0.7` |
| `--fivemc-code` | Code used for 5mC in modkit files. | `m` |

Calls with a score between the two thresholds are ignored. Each probe is assigned methylated or unmethylated by majority vote of the remaining calls within the margin.

### `single` and `multi`

| Option | Description | Default |
|---|---|---|
| `-i`, `--input_file` (`single`)<br>`-i`, `--input_dir` (`multi`) | Path to the input file (`single`) or to a directory of input files (`multi`). | — (required) |
| `-m`, `--model` | Path to the model zip file, e.g. `cns-v2.4b.zip`. | — (required) |
| `-o`, `--output_file` | Path of the output file (CSV for `single`, tab-separated for `multi`). | — (required) |
| `-f`, `--file_type` | Input file type: `bed` or `csv`. In `multi` mode, only files with the matching extension are used. | — (required) |
| `-c`, `--num_CpGs` | Number of CpG sites to randomly sample. | All available sites |
| `-s`, `--num_samplings` | Number of random samplings of size `num_CpGs`. | `1` |

### Input formats

**bed**: the output of `sturgeon-v2 inputtobed` (or of the original `sturgeon inputtobed`). A whitespace-separated file with a header; Sturgeon V2 uses the `probe_id` and `methylation_call` columns (`1` = methylated, `0` = unmethylated).

```
chrom	chromStart	chromEnd	methylation_call	probe_id
1	380981	380982	1	cg01394461
1	381818	381819	1	cg09555360
```

**csv**: a comma-separated file **without a header** and with two columns: probe ID and methylation call (`1` = methylated, `0` = unmethylated).

```
cg01394461,1
cg09555360,1
cg14057946,0
```

In both formats, probes that are not listed are treated as not measured.

## Output

### `single`

A comma-separated file with one row per sampling (so a single row when `-s` is not set), and the following columns:

| Column | Description |
|---|---|
| one column per class | The calibrated score for that class. Scores in each row sum to 1. The class with the highest score is the prediction. |
| `probes` | Number of CpG sites used for the prediction. |

Scores are calibrated for sequencing depth: the model output is scaled with a temperature that depends on the number of CpG sites used, before being converted to probabilities.

<!-- TODO: add a short example table with real class names from cns-v2.4b -->

### `multi`

A tab-separated file with one row per input file:

| Column | Description |
|---|---|
| `sample` | Input file name without its extension. |
| `class` | The class with the highest score, averaged over all samplings. |
| `score` | The averaged score of that class. |
| `purity` | Reserved; currently left empty. |


