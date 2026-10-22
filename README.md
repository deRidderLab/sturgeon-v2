# Sturgeon V2

**Intraoperative CNS tumour classification from nanopore methylation data.**

<img width="1362" height="426" alt="image" src="https://github.com/user-attachments/assets/f90597f1-641d-4cca-beb8-ccbc02d8b4ff" />


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
- [modkit](https://github.com/nanoporetech/modkit), for extracting methylation calls from BAM files

The original Sturgeon package is **not** required: preprocessing (`inputtobed`) is included in Sturgeon V2.

### Install Sturgeon V2

Install directly from GitHub over SSH (requires an SSH key linked to your GitHub account):

```bash
pip install git+ssh://github.com/deRidderLab/sturgeon-v2.git
```

Alternatively, install over HTTPS:

```bash
pip install git+https://github.com/deRidderLab/sturgeon-v2.git
```

Or install from a local clone:

```bash
git clone https://github.com/deRidderLab/sturgeon-v2.git
cd sturgeon-v2
pip install .
```

Verify the installation:

```bash
sturgeon-v2 --help
```

## Models

Models are distributed as zip files and are passed to `sturgeon-v2` with the `-m/--model` option, both for preprocessing (`inputtobed`) and for prediction (`single`/`multi`). Each model zip contains the probe locations (`probe.bed`) that the model was trained on, so preprocessing and prediction always use the same set of probes. Download a model and store it in a location of your choice; it does not need to be unpacked.

| Model | Status | Download |
|---|---|---|
| `cns-v2` | **Current, recommended** | [zip file](https://zenodo.org/records/22899108) |
| `cns-v2.4a` | Deprecated | [zip file](https://drive.google.com/file/d/1opEUbGWtUCONWhlfrv9NfSuGM8pXj-D8/view?usp=share_link) |

<!-- TODO: describe the classification scheme of cns-v2.4b (number of classes, control classes, merged classes) and what changed between v2.4a and v2.4b -->

> **Use the same model for preprocessing and prediction.** Pass the same model zip to `inputtobed` and to `single`/`multi`.

## Workflow overview

```
Nanopore BAM (with 5mC calls, aligned to T2T-CHM13v2.0)
        │
        │  modkit adjust-mods / modkit extract        (Step 1a)
        ▼
modkit extract output (.txt)
        │
        │  sturgeon-v2 inputtobed -m MODEL            (Step 1b)
        ▼
Methylation calls on array probe sites (.bed)
        │
        │  sturgeon-v2 single / multi -m MODEL        (Step 2)
        ▼
Class scores per sample
```

## Step 1: Preprocessing

Preprocessing is done with `sturgeon-v2 inputtobed`. This command is ported from the `inputtobed` utility of the original [Sturgeon](https://github.com/UMCUGenetics/sturgeon) and uses the same mapping logic and default settings, so `.bed` files previously created with `sturgeon inputtobed` can still be used with Sturgeon V2. The probe locations are read from the model zip (`-m`).

> ⚠️ **Reference genome: only T2T-CHM13v2.0 is supported.**
> Sturgeon V2 does **not** support `hg38`. Reads must be aligned to the [T2T-CHM13v2.0 reference genome](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz), which the probe locations in the model zip are based on. If your data is aligned to `hg38`, re-align it to T2T-CHM13v2.0 before preprocessing. Incorrectly aligned data will give wrong results without raising an error.

### 1a. Extract methylation calls with modkit

Starting from a BAM file with modified base calls (from Dorado or Guppy), aligned to T2T-CHM13v2.0:

**(Optional) Combine 5hmC with 5mC.** If your BAM contains both 5mC and 5hmC calls, merge them into a single 5mC score. Both modifications look alike on bisulfite-based methylation arrays, which the reference data was generated with, so combining them best matches the training data.

```bash
modkit adjust-mods --convert h m input.bam adjusted.bam
```

**Extract per-read calls:**

```bash
mkdir -p modkit
modkit extract full adjusted.bam modkit/sample1.txt
```

Older versions of modkit use `modkit extract` without `full`. Consult the [modkit documentation](https://nanoporetech.github.io/modkit/quick_start.html) if the command has changed.

Name the output file after the sample: `inputtobed` uses the file name for the bed file it writes, and in `multi` mode that name becomes the sample name.

### 1b. Convert to bed with `sturgeon-v2 inputtobed`

`inputtobed` takes `modkit extract` output and the model zip. It can convert one sample or a whole directory of samples.

**One sample** (for `single`): pass the modkit file; `-o` is the bed file to write.

```bash
sturgeon-v2 inputtobed -i modkit/sample1.txt -o bed_files/sample1.bed -m models/cns-v2.4b.zip
```

**Several samples** (for `multi`): pass a directory; every `.txt` (or `.tsv`) modkit file in it is converted separately to a bed file with the same name, so `sample1.txt` becomes `sample1.bed` in the output directory. Files are never merged.

```bash
sturgeon-v2 inputtobed -i modkit/ -o bed_files/ -m models/cns-v2.4b.zip
```

Only the bed files are written; no intermediate files are created. Re-running overwrites existing bed files with the same name.

> **One file per sample.** In directory mode, each modkit file is treated as a separate sample. If the modkit output of one sample is split over several files, combine them into one file first (keeping a single header line), for example:
> ```bash
> head -n 1 part1.txt > sample1.txt
> tail -q -n +2 part*.txt >> sample1.txt
> ```

## Step 2: Prediction

Sturgeon V2 has two prediction modes:

- `single` predicts on one `.bed` or `.csv` file.
- `multi` predicts on every `.bed` or `.csv` file in a directory, treating each file as an independent sample.

### Predict a single sample

```bash
sturgeon-v2 single \
    -i bed_files/sample1.bed \
    -m models/cns-v2.4b.zip \
    -f bed \
    -o results/sample1 \
    -p
```

This writes `results/sample1/result.txt` and, because of `-p`, `results/sample1/plot.png`. The predicted class, its score and the number of CpG sites used are also printed in the terminal.

### Predict a directory of samples

Use the output directory of `inputtobed` in directory mode. The file name (without extension) is used as the sample name in the output.

```bash
sturgeon-v2 multi \
    -i bed_files/ \
    -m models/cns-v2.4b.zip \
    -f bed \
    -o results/all_samples \
    -p
```

This writes `results/all_samples/result.txt`, one file covering all samples, and, because of `-p`, one plot per sample in `results/all_samples/plots/` (`sample1.png`, `sample2.png`, …).

### Plots

With the `-p` flag, `single` and `multi` also save plots in the output directory: one per sample, with a vertical bar for every class in the order of the model's `classification_system.yaml`. Each bar is the score of that class, averaged over all samplings; the title shows the predicted class, its score and the number of CpG sites used. Class colors are taken from the same yaml file; classes without a color are shown in grey.

- `single`: `<output_dir>/plot.png`.
- `multi`: `<output_dir>/plots/<sample>.png`, one per input file.

### Simulating lower coverage (subsampling)

The `-c` and `-s` options of `single` and `multi` randomly subsample the available CpG sites. This is useful for benchmarking how the classifier behaves at lower sequencing depth, for example to mimic an early timepoint during surgery. Sampling is seeded, so results are reproducible.

The following example draws 10 random subsets of 1,000 CpG sites and predicts on each:

```bash
sturgeon-v2 single \
    -i bed_files/sample1.bed \
    -m models/cns-v2.4b.zip \
    -f bed \
    -c 1000 \
    -s 10 \
    -o results/sample1_1000cpgs \
    -p
```

If `-c` is larger than the number of available CpG sites, all sites are used. For routine clinical use, leave `-c` unset so that all available CpG sites are used.

## Command-line reference

`sturgeon-v2` has three subcommands:

```
sturgeon-v2 inputtobed -i INPUT_PATH -o OUTPUT_PATH -m MODEL
sturgeon-v2 single     -i INPUT_FILE -m MODEL -o OUTPUT_DIR -f FILE_TYPE [-c NUM_CPGS] [-s NUM_SAMPLINGS] [-p]
sturgeon-v2 multi      -i INPUT_DIR  -m MODEL -o OUTPUT_DIR -f FILE_TYPE [-c NUM_CPGS] [-s NUM_SAMPLINGS] [-p]
```

Run `sturgeon-v2 <subcommand> --help` for the full list of options.

### `inputtobed`

| Option | Description | Default |
|---|---|---|
| `-i`, `--input-path` | A `modkit extract` file, or a directory of them (`.txt` or `.tsv`). | — (required) |
| `-o`, `--output-path` | The output `.bed` file when `-i` is a file; the output directory when `-i` is a directory (one `<name>.bed` per input file). | — (required) |
| `-m`, `--model` | Path to the model zip file, e.g. `cns-v2.4b.zip`. The probe locations (`probe.bed`) are read from it. Use the same model as for prediction. | — (required) |

The conversion uses fixed settings, the defaults of the original Sturgeon: for each probe, all 5mC calls (modkit code `m`) within 25 bp are considered. Calls with a score above 0.7 count as methylated, calls below 0.3 as unmethylated, and calls in between are ignored. Each probe is assigned methylated or unmethylated by majority vote; probes with a tie or without calls are left out.

### `single` and `multi`

| Option | Description | Default |
|---|---|---|
| `-i`, `--input_file` (`single`)<br>`-i`, `--input_dir` (`multi`) | Path to the input file (`single`) or to a directory of input files (`multi`). | — (required) |
| `-m`, `--model` | Path to the model zip file, e.g. `cns-v2.4b.zip`. | — (required) |
| `-o`, `--output_dir` | Output directory, created if needed. The results are written to `result.txt` in it; see [Output](#output). | — (required) |
| `-f`, `--file_type` | Input file type: `bed` or `csv`. In `multi` mode, only files with the matching extension are used. | — (required) |
| `-c`, `--num_CpGs` | Number of CpG sites to randomly sample. | All available sites |
| `-s`, `--num_samplings` | Number of random samplings of size `num_CpGs`. | `1` |
| `-p`, `--plot` | Flag: also save plots in the output directory, `plot.png` (`single`) or `plots/<sample>.png` (`multi`); see [Plots](#plots). | No plots |

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

`single` and `multi` write their results to `result.txt` in the output directory (`-o`). With `-p`, the plots are written to the same directory:

```
OUTPUT_DIR/                 OUTPUT_DIR/
├── result.txt              ├── result.txt
└── plot.png   (single)     └── plots/       (multi)
                                ├── sample1.png
                                └── sample2.png
```

Existing files with the same name are overwritten, so use a new output directory per run if you want to keep earlier results.

Scores are calibrated for sequencing depth: the model output is scaled with a temperature that depends on the number of CpG sites used, before being converted to probabilities. The scores of all classes sum to 1.

### `single`

`result.txt` is a tab-separated file with the raw scores: one row per sampling (so a single row when `-s` is not set), and the following columns:

| Column | Description |
|---|---|
| one column per class | The calibrated score for that class, in the order of the model's `classification_system.yaml`. |
| `probes` | Number of CpG sites used for the prediction. |

The predicted class (the class with the highest score, averaged over all samplings), its score and the number of CpG sites used are also printed in the terminal.

### `multi`

`result.txt` is a tab-separated file with one summary row per input file:

| Column | Description |
|---|---|
| `sample` | Input file name without its extension. |
| `class` | The predicted class: the class with the highest score, averaged over all samplings. |
| `score` | The averaged score of that class. |

Example:

```
sample	class	score
sample1	MB_SHH	0.9731
sample2	GBM_RTK2	0.8812
```
