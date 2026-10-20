import argparse
import io
import os
import re
import shutil
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import polars as pl
import yaml
import onnxruntime as ort
from tqdm import tqdm

from sturgeon_v2.plotting import get_class_colors, plot_scores


SUMMARY_COLUMNS = ["sample", "class", "score"]


# ---------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------
def softmax(x, axis=1):
    x = x - np.max(x, axis=axis, keepdims=True)
    exp_x = np.exp(x)
    return exp_x / np.sum(exp_x, axis=axis, keepdims=True)


# ---------------------------------------------------------------------
# Model loading
# ---------------------------------------------------------------------
@dataclass
class Model:
    sess: ort.InferenceSession
    input_name: str
    output_name: str
    class_names: List[str]
    colors: Dict[str, str]
    probe_df: pl.DataFrame
    bin_centers: np.ndarray
    temps: np.ndarray


def load_model(model_zip_path: str, tmpdir: str) -> Model:
    with zipfile.ZipFile(model_zip_path, "r") as zf:
        zf.extractall(tmpdir)

    model_root = os.path.join(tmpdir, '.'.join(os.path.basename(model_zip_path).split('.')[:-1]))
    model_path = os.path.join(model_root, "model.onnx")

    if not os.path.exists(model_path):
        raise FileNotFoundError(f"model.onnx not found at {model_path}")

    with open(os.path.join(model_root, "classification_system.yaml"), "r") as f:
        classification_file = yaml.safe_load(f)

    probe_df = pl.read_csv(
        os.path.join(model_root, "probe.bed"),
        separator="\t",
    )

    print("-------------------------------")
    print("Loading ONNX model")
    print("-------------------------------")

    so = ort.SessionOptions()
    so.intra_op_num_threads = 1
    so.inter_op_num_threads = 1

    sess = ort.InferenceSession(
        model_path,
        sess_options=so,
        providers=["CPUExecutionProvider"],
    )

    return Model(
        sess=sess,
        input_name=sess.get_inputs()[0].name,
        output_name=sess.get_outputs()[0].name,
        class_names=list(classification_file["encoder"]["type"].keys()),
        colors=get_class_colors(classification_file),
        probe_df=probe_df,
        bin_centers=np.load(os.path.join(model_root, "bin_centers.npy")),
        temps=np.load(os.path.join(model_root, "temps.npy")),
    )


# ---------------------------------------------------------------------
# Input reading
# ---------------------------------------------------------------------
def read_input(input_file: str, file_type: str, probe_df: pl.DataFrame) -> np.ndarray:
    """Read a bed or csv file and return the model input vector
    (1 = methylated, -1 = unmethylated, 0 = not measured)."""
    if file_type == "csv":
        bed_file = pl.read_csv(
            input_file,
            separator=",",
            has_header=False,
            new_columns=["probe_id", "methylation_call"],
        )
    elif file_type == "bed":
        with open(input_file, "r") as f:
            raw_text = f.read()
        # collapse any run of whitespace (spaces, tabs, mixed) into a single tab
        normalized = re.sub(r"[ \t]+", "\t", raw_text)
        bed_file = pl.read_csv(io.StringIO(normalized), separator="\t")
    else:
        raise ValueError(f"Unsupported file type: {file_type}")

    bed_file = bed_file.select(
        [
            pl.col("probe_id").cast(pl.Utf8),
            pl.col("methylation_call")
            .cast(pl.Int32, strict=False)
            .alias("methylation_call"),
        ]
    )

    bed_file = bed_file.with_columns(
        methylation_call=pl.when(pl.col("methylation_call") == 0)
        .then(-1)
        .otherwise(pl.col("methylation_call"))
    )

    nn_input = probe_df.join(
        bed_file,
        left_on="name",
        right_on="probe_id",
        how="left",
        validate="1:1",
    ).with_columns(pl.all().fill_null(0))

    return nn_input["methylation_call"].to_numpy().astype(np.float32)


# ---------------------------------------------------------------------
# Prediction
# ---------------------------------------------------------------------
def predict_scores(
    model: Model,
    nn_input: np.ndarray,
    n: Optional[str],
    num_samples: int,
) -> Tuple[np.ndarray, int, bool]:
    """Run the (sub)sampled predictions for one sample.

    Returns the calibrated scores (num_samples x num_classes), the number of
    CpG sites used, and whether all available sites were used.
    """
    input_size = len(nn_input)
    positions = np.nonzero(nn_input)[0]

    use_all = n is None or int(n) > len(positions)
    n_used = len(positions) if use_all else int(n)

    result = np.zeros((num_samples, len(model.class_names)))
    for i in range(num_samples):
        rng = np.random.default_rng(seed=i)
        selected = rng.choice(positions, size=n_used, replace=False)

        new_input = np.zeros_like(nn_input)
        new_input[selected] = nn_input[selected]
        new_input = new_input.reshape(1, input_size)

        outputs = model.sess.run(
            [model.output_name],
            {model.input_name: new_input},
        )[0]

        index = np.argmin(np.abs(model.bin_centers - n_used))
        temperature = model.temps[index]
        outputs = outputs / np.exp(temperature)

        outputs = softmax(outputs, axis=1)
        result[i, :] = outputs[0]

    return result, n_used, use_all


def summarise(sample_name: str, class_names: List[str], result: np.ndarray) -> dict:
    """Average the scores over the samplings and keep the best class."""
    avg_probs = np.mean(result, axis=0)
    best_idx = int(np.argmax(avg_probs))
    return {
        "sample": sample_name,
        "class": class_names[best_idx],
        "score": float(avg_probs[best_idx]),
    }


def _ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)


# ---------------------------------------------------------------------
# Single File Prediction
# ---------------------------------------------------------------------
def predict(
    input_file,
    model_zip_path,
    n,
    output_file,
    num_samples,
    file_type,
    plot_file=None,
):
    if num_samples is None:
        num_samples = 1

    tmpdir = tempfile.mkdtemp(prefix="sturgeon_onnx_")

    try:
        model = load_model(model_zip_path, tmpdir)
        nn_input = read_input(input_file, file_type, model.probe_df)
        sample_name = Path(input_file).stem

        print("-------------------------------")
        print("Running subsample predictions")
        print("-------------------------------")

        result, n_used, use_all = predict_scores(model, nn_input, n, num_samples)
        if use_all:
            print("Using all positions")

        summary = summarise(sample_name, model.class_names, result)

        # raw scores: one row per sampling, one column per class
        _ensure_parent_dir(output_file)
        df = pd.DataFrame(result, columns=model.class_names)
        df["probes"] = n_used
        df.to_csv(output_file, index=False)

        print("-------------------------------")
        print(f"Sample:          {summary['sample']}")
        print(f"Predicted class: {summary['class']}")
        print(f"Score:           {summary['score']:.4f}"
              + (f" (mean over {num_samples} samplings)" if num_samples > 1 else ""))
        print(f"CpG sites used:  {n_used}")
        print("-------------------------------")
        print(f"Results saved to {output_file}")

        if plot_file:
            # plots are always stored as .png
            if Path(plot_file).suffix.lower() != ".png":
                plot_file = str(Path(plot_file).with_suffix(".png"))
            _ensure_parent_dir(plot_file)
            plot_scores(
                sample_name=sample_name,
                class_names=model.class_names,
                scores=result.mean(axis=0),
                colors=model.colors,
                n_cpgs=n_used,
                output_file=plot_file,
            )
            print(f"Plot saved to {plot_file}")

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ---------------------------------------------------------------------
# Directory Multi-Prediction
# ---------------------------------------------------------------------
def predict_directory(
    input_dir,
    model_zip_path,
    n,
    output_file,
    num_samples,
    file_type,
    plot_dir=None,
):
    if num_samples is None:
        num_samples = 1

    tmpdir = tempfile.mkdtemp(prefix="sturgeon_onnx_")

    try:
        model = load_model(model_zip_path, tmpdir)
        if plot_dir:
            os.makedirs(plot_dir, exist_ok=True)

        print("-------------------------------")
        print(f"Running predictions on directory: {input_dir}")
        print("-------------------------------")

        valid_exts = {"csv": ".csv", "bed": ".bed"}
        if file_type not in valid_exts:
            raise ValueError(f"Unsupported file type: {file_type}")
        ext = valid_exts[file_type]

        results_list = []

        for filename in tqdm(sorted(os.listdir(input_dir))):
            if not filename.endswith(ext):
                continue

            input_file = os.path.join(input_dir, filename)
            sample_name = filename[:-len(ext)]

            nn_input = read_input(input_file, file_type, model.probe_df)
            result, n_used, _ = predict_scores(model, nn_input, n, num_samples)
            summary = summarise(sample_name, model.class_names, result)

            results_list.append(summary)

            if plot_dir:
                plot_scores(
                    sample_name=sample_name,
                    class_names=model.class_names,
                    scores=result.mean(axis=0),
                    colors=model.colors,
                    n_cpgs=n_used,
                    output_file=os.path.join(plot_dir, f"{sample_name}.png"),
                )
            print(f"Processed {filename} -> {summary['class']} ({summary['score']:.4f})")

        if not results_list:
            raise ValueError(f"No {ext} files found in {input_dir}")

        _ensure_parent_dir(output_file)
        df_out = pd.DataFrame(results_list, columns=SUMMARY_COLUMNS)
        df_out.to_csv(output_file, sep="\t", index=False)
        print(f"\nAll multi-results saved to {output_file}")
        if plot_dir:
            print(f"Plots saved to {plot_dir}")

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ---------------------------------------------------------------------
# CLI & Entrypoint
# ---------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Sturgeon v2: Intraoperative methylation classification")
    
    # Create subparsers for 'single' and 'multi' modes
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Parser for `sturgeon-v2 single`
    parser_single = subparsers.add_parser("single", help="Predict on a single .bed or .csv file")
    parser_single.add_argument("-i", "--input_file", required=True, help="Input file path")
    parser_single.add_argument("-m", "--model", required=True, help="Path to the model zip")
    parser_single.add_argument("-c", "--num_CpGs", help="Number of CpGs to use")
    parser_single.add_argument("-o", "--output_file", required=True, help="Path to save the output csv (scores of all classes per sampling)")
    parser_single.add_argument("-s", "--num_samplings", help="Number of subsamplings")
    parser_single.add_argument("-f", "--file_type", required=True, choices=["bed", "csv"], help="Type of file")
    parser_single.add_argument("-p", "--plot", default=None,
                               help="Save a plot of the scores of all classes to this .png file")

    # Parser for `sturgeon-v2 multi`
    parser_multi = subparsers.add_parser("multi", help="Predict on a directory of .bed or .csv files")
    parser_multi.add_argument("-i", "--input_dir", required=True, help="Input directory path")
    parser_multi.add_argument("-m", "--model", required=True, help="Path to the model zip")
    parser_multi.add_argument("-c", "--num_CpGs", help="Number of CpGs to use")
    parser_multi.add_argument("-o", "--output_file", required=True, help="Path to save the combined output tsv")
    parser_multi.add_argument("-s", "--num_samplings", help="Number of subsamplings")
    parser_multi.add_argument("-f", "--file_type", required=True, choices=["bed", "csv"], help="Type of files")
    parser_multi.add_argument("-p", "--plot_dir", default=None,
                              help="Directory to save one plot per sample (<sample>.png) with the scores of all classes")

    # Parser for `sturgeon-v2 inputtobed`
    parser_bed = subparsers.add_parser(
        "inputtobed",
        help="Convert modkit extract output to bed file(s) for single or multi",
    )
    parser_bed.add_argument("-i", "--input-path", required=True,
                            help="A modkit extract file (gives one bed, for 'single'), or a directory of "
                                 "modkit extract files (gives one bed per file, for 'multi')")
    parser_bed.add_argument("-o", "--output-path", required=True,
                            help="Output bed file when the input is a file; output directory when the "
                                 "input is a directory")
    parser_bed.add_argument("-m", "--model", required=True,
                            help="Path to the model zip; its probe.bed is used for the probe locations "
                                 "(use the same model as for prediction)")

    args = parser.parse_args()

    if args.command == "inputtobed":
        import logging
        from sturgeon_v2.inputtobed import modkit_to_bed

        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
        print("---------------------------------------------------------------------")
        print("Converting modkit methylation calls to bed")
        print("---------------------------------------------------------------------")
        bed_files = modkit_to_bed(
            input_path=args.input_path,
            output_path=args.output_path,
            model_zip=args.model,
        )
        if os.path.isdir(args.input_path):
            print(f"{len(bed_files)} bed files saved to {args.output_path}")
            print(f"Run prediction with: sturgeon-v2 multi -f bed -i {args.output_path} "
                  f"-m {args.model} -o OUTPUT_TSV")
        else:
            print(f"Bed file saved to {bed_files[0]}")
            print(f"Run prediction with: sturgeon-v2 single -f bed -i {bed_files[0]} "
                  f"-m {args.model} -o OUTPUT_TSV")
        return

    num_samples = int(args.num_samplings) if args.num_samplings else 1

    print("---------------------------------------------------------------------")
    print(f"Running intraoperative methylation based classification (ONNX) - Mode: {args.command.upper()}")
    print("---------------------------------------------------------------------")

    if args.command == "single":
        predict(
            input_file=args.input_file,
            model_zip_path=args.model,
            n=args.num_CpGs,
            output_file=args.output_file,
            num_samples=num_samples,
            file_type=args.file_type,
            plot_file=args.plot,
        )
    elif args.command == "multi":
        predict_directory(
            input_dir=args.input_dir,
            model_zip_path=args.model,
            n=args.num_CpGs,
            output_file=args.output_file,
            num_samples=num_samples,
            file_type=args.file_type,
            plot_dir=args.plot_dir,
        )

def entrypoint():
    try:
        main()
    except Exception:
        import traceback
        print("\n❌ ERROR: An exception occurred in sturgeon-v2 CLI\n")
        traceback.print_exc()
        exit(1)

if __name__ == "__main__":
    entrypoint()
