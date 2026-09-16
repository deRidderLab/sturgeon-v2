# Sturgeon V2

**Intraoperative CNS tumour classification from nanopore methylation data.**

Sturgeon V2 is the successor of [Sturgeon](https://github.com/UMCUGenetics/sturgeon), a neural network classifier for central nervous system (CNS) tumours based on sparse DNA methylation profiles, such as those obtained after only a few minutes of nanopore sequencing during surgery. This repository provides a pip-installable package to run Sturgeon V2 in a clinical setting.

## Table of contents
- [Installation](#installation)
- [Models](#models)
- [Workflow overview](#workflow-overview)
- [Step 1: Preprocessing (using Sturgeon)](#step-1-preprocessing-using-sturgeon)
- [Step 2: Prediction (using Sturgeon V2)](#step-2-prediction-using-sturgeon-v2)
- [Command-line reference](#command-line-reference)
- [Output](#output)
- [Interpreting results](#interpreting-results)

## Installation

### Requirements

- Python <!-- TODO: e.g. >= 3.9 -->
- Linux or macOS <!-- TODO: confirm supported platforms -->
- The original [Sturgeon](https://github.com/UMCUGenetics/sturgeon) package, for preprocessing (see [Step 1](#step-1-preprocessing-using-sturgeon))
- [modkit](https://github.com/nanoporetech/modkit) (recommended), for extracting methylation calls from BAM files

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

### Install Sturgeon (for preprocessing)

Preprocessing is done with the original Sturgeon package. Follow the installation instructions in the [Sturgeon repository](https://github.com/UMCUGenetics/sturgeon#installation). You may install it in the same environment as Sturgeon V2 or in a separate one.

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
        │  sturgeon inputtobed                        (Step 1b)
        ▼
Methylation calls on array probe sites (.bed)
        │
        │  sturgeon-v2 single / multi                 (Step 2)
        ▼
Class scores per sample
```

## Step 1: Preprocessing (using Sturgeon)

Sturgeon V2 uses exactly the same preprocessing as the original Sturgeon. The steps are summarised below; for full details see the [Sturgeon README](https://github.com/UMCUGenetics/sturgeon#preparing-data-into-the-right-format-inputtobed).

> ⚠️ **Reference genome: only T2T-CHM13v2.0 is supported.**
> Sturgeon V2 does **not** support `hg38`. Reads must be aligned to the [T2T-CHM13v2.0 reference genome](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz). Do **not** pass `--reference-genome hg38` to `sturgeon inputtobed`. If your data is aligned to `hg38`, re-align it to T2T-CHM13v2.0 before preprocessing. Incorrectly aligned data will give wrong results without raising an error.

### 1a. Extract methylation calls with modkit

Starting from a BAM file with modified base calls (from Dorado or Guppy), aligned to T2T-CHM13v2.0:

**(Optional) Combine 5hmC with 5mC.** If your BAM contains both 5mC and 5hmC calls, merge them into a single 5mC score. Both modifications look alike on bisulfite-based methylation arrays, which the reference data was generated with, so combining them best matches the training data.

```bash
modkit adjust-mods --convert h m input.bam adjusted.bam
```

**Extract per-read calls:**

```bash
modkit extract full adjusted.bam modkit_output.txt
```

Older versions of modkit use `modkit extract` without `full`. Consult the [modkit documentation](https://nanoporetech.github.io/modkit/quick_start.html) if the command has changed.

Alternatively, `modkit pileup` can be used to produce a single file of aggregated calls per genomic position.

### 1b. Convert to bed with `sturgeon inputtobed`

For `modkit extract` output:

```bash
sturgeon inputtobed -i MODKIT_OUTPUT_DIR -o OUTPUT_DIR -s modkit
```

For `modkit pileup` output:

```bash
sturgeon inputtobed -i MODKIT_PILEUP_FILE -o OUTPUT_FILE -s modkit_pileup
```

The resulting `.bed` file (e.g. `merged_probes_methyl_calls.bed`) is the input for Sturgeon V2.

> **One sample per run.** `inputtobed` merges all methylation calls from all files in the input directory into a single output. If you have multiple samples, process each sample separately into its own output directory.

## Step 2: Prediction (using Sturgeon V2)

Sturgeon V2 has two modes:

- `single` predicts on one `.bed` or `.csv` file.
- `multi` predicts on every `.bed` or `.csv` file in a directory, treating each file as an independent sample.

### Predict a single sample

```bash
sturgeon-v2 single \
    -i sample1/merged_probes_methyl_calls.bed \
    -m models/cns-v2.4b.zip \
    -f bed \
    -o results/sample1_prediction.csv
```

### Predict a directory of samples

```bash
sturgeon-v2 multi \
    -i bed_files/ \
    -m models/cns-v2.4b.zip \
    -f bed \
    -o results/
```

<!-- TODO: confirm whether -o in multi mode expects a directory or a single file -->

### Simulating lower coverage (subsampling)

The `-c` and `-s` options randomly subsample the available CpG sites. This is useful for benchmarking how the classifier behaves at lower sequencing depth, for example to mimic an early timepoint during surgery.

The following example draws 10 random subsets of 1,000 CpG sites and predicts on each:

```bash
sturgeon-v2 single \
    -i sample1/merged_probes_methyl_calls.bed \
    -m models/cns-v2.4b.zip \
    -f bed \
    -c 1000 \
    -s 10 \
    -o results/sample1_1000cpgs.csv
```

For routine clinical use, leave `-c` unset so that all available CpG sites are used.

## Command-line reference

```
sturgeon-v2 {single,multi} -i INPUT -m MODEL -o OUTPUT [-f FILE_TYPE] [-c NUM_CPGS] [-s NUM_SAMPLINGS]
```

| Option | Description | Default |
|---|---|---|
| `single` / `multi` | Predict on a single file, or on all files in a directory. | — |
| `-i`, `--input_file` | Path to the input file (`single`) or directory (`multi`). | — |
| `-m`, `--model` | Path to the model zip file, e.g. `cns-v2.4b.zip`. | — |
| `-o`, `--output_file` | Path to write the predictions to. | — |
| `-f`, `--file_type` | Input file type: `bed` or `csv`. | <!-- TODO --> |
| `-c`, `--num_CpGs` | Number of CpG sites to randomly sample. | All available sites |
| `-s`, `--num_samplings` | Number of random samplings of size `num_CpGs`. | `1` |

### Input formats

**bed**: the output of `sturgeon inputtobed` (see [Step 1](#step-1-preprocessing-using-sturgeon)).

**csv**: <!-- TODO: describe required columns, e.g. probe ID and methylation call, with a small example -->

## Output

<!-- TODO: describe the output file: columns (one per class?), the calibrated score, which column holds the top prediction, number of CpGs used, and one row per sampling when -s > 1. A short example table would help. -->
