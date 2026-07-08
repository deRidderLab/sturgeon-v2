import argparse
import io
import zipfile
import tempfile
import os
import shutil

import numpy as np
import pandas as pd
import polars as pl
import yaml
import onnxruntime as ort
from tqdm import tqdm


# ---------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------
def softmax(x, axis=1):
    x = x - np.max(x, axis=axis, keepdims=True)
    exp_x = np.exp(x)
    return exp_x / np.sum(exp_x, axis=axis, keepdims=True)


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
):
    if num_samples is None:
        num_samples = 1

    tmpdir = tempfile.mkdtemp(prefix="sturgeon_onnx_")

    try:
        with zipfile.ZipFile(model_zip_path, "r") as zf:
            zf.extractall(tmpdir)

        model_root = os.path.join(tmpdir, '.'.join(os.path.basename(model_zip_path).split('.')[:-1]))
        model_path = os.path.join(model_root, "model.onnx")

        if not os.path.exists(model_path):
            raise FileNotFoundError(f"model.onnx not found at {model_path}")

        with open(os.path.join(model_root, "classification_system.yaml"), "r") as f:
            classification_file = yaml.safe_load(f)

        class_names = list(classification_file["encoder"]["type"].keys())
        output_size = len(class_names)

        probe_df = pl.read_csv(
            os.path.join(model_root, "probe.bed"),
            separator="\t",
        )
        input_size = len(probe_df)

        if file_type == "csv":
            bed_file = pl.read_csv(
                input_file,
                separator=",",
                has_header=False,
                new_columns=["probe_id", "methylation_call"],
            )
        elif file_type == "bed":
            bed_file = pl.read_csv(input_file, separator="\t")
        elif file_type == "bed_carlo":
            bed_file = pl.read_csv(input_file, separator=" ")
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

        nn_input = nn_input["methylation_call"].to_numpy().astype(np.float32)

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

        input_name = sess.get_inputs()[0].name
        output_name = sess.get_outputs()[0].name

        bin_centers = np.load(os.path.join(model_root, "bin_centers.npy"))
        temps = np.load(os.path.join(model_root, "temps.npy"))

        print("-------------------------------")
        print("Running subsample predictions")
        print("-------------------------------")

        result = np.zeros((num_samples, output_size))
        positions = np.nonzero(nn_input)[0]

        for i in range(num_samples):
            if n is None or int(n) > len(positions):
                if i == 0:
                    print("-------------------------------")
                    print("Using all positions")
                    print("-------------------------------")
                n_used = len(positions)
            else:
                n_used = int(n)

            rng = np.random.default_rng(seed=i)
            selected = rng.choice(positions, size=n_used, replace=False)

            new_input = np.zeros_like(nn_input)
            new_input[selected] = nn_input[selected]
            new_input = new_input.reshape(1, input_size)

            outputs = sess.run(
                [output_name],
                {input_name: new_input},
            )[0]

            index = np.argmin(np.abs(bin_centers - n_used))
            temperature = temps[index]
            outputs = outputs / np.exp(temperature)

            outputs = softmax(outputs, axis=1)
            result[i, :] = outputs[0]

        df = pd.DataFrame(result, columns=class_names)
        df["probes"] = n_used
        df.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")

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
):
    if num_samples is None:
        num_samples = 1

    tmpdir = tempfile.mkdtemp(prefix="sturgeon_onnx_")

    try:
        with zipfile.ZipFile(model_zip_path, "r") as zf:
            zf.extractall(tmpdir)

        model_root = os.path.join(tmpdir, '.'.join(os.path.basename(model_zip_path).split('.')[:-1]))
        model_path = os.path.join(model_root, "model.onnx")

        if not os.path.exists(model_path):
            raise FileNotFoundError(f"model.onnx not found at {model_path}")

        with open(os.path.join(model_root, "classification_system.yaml"), "r") as f:
            classification_file = yaml.safe_load(f)

        class_names = list(classification_file["encoder"]["type"].keys())
        output_size = len(class_names)

        probe_df = pl.read_csv(
            os.path.join(model_root, "probe.bed"),
            separator="\t",
        )
        input_size = len(probe_df)

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

        input_name = sess.get_inputs()[0].name
        output_name = sess.get_outputs()[0].name

        bin_centers = np.load(os.path.join(model_root, "bin_centers.npy"))
        temps = np.load(os.path.join(model_root, "temps.npy"))

        print("-------------------------------")
        print(f"Running predictions on directory: {input_dir}")
        print("-------------------------------")

        valid_exts = {"csv": ".csv", "bed": ".bed", "bed_carlo": ".bed"}
        ext = valid_exts.get(file_type, ".bed")
        
        results_list = []

        for filename in tqdm(sorted(os.listdir(input_dir))):
            if not filename.endswith(ext):
                continue
                
            file_path = os.path.join(input_dir, filename)
            sample_name = filename[:-len(ext)]

            if file_type == "csv":
                bed_file = pl.read_csv(
                    file_path,
                    separator=",",
                    has_header=False,
                    new_columns=["probe_id", "methylation_call"],
                )
            elif file_type == "bed":
                bed_file = pl.read_csv(file_path, separator="\t")
            elif file_type == "bed_carlo":
                bed_file = pl.read_csv(file_path, separator=" ")

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

            nn_input = nn_input["methylation_call"].to_numpy().astype(np.float32)

            result = np.zeros((num_samples, output_size))
            positions = np.nonzero(nn_input)[0]

            for i in range(num_samples):
                if n is None or int(n) > len(positions):
                    n_used = len(positions)
                else:
                    n_used = int(n)

                rng = np.random.default_rng(seed=i)
                selected = rng.choice(positions, size=n_used, replace=False)

                new_input = np.zeros_like(nn_input)
                new_input[selected] = nn_input[selected]
                new_input = new_input.reshape(1, input_size)

                outputs = sess.run(
                    [output_name],
                    {input_name: new_input},
                )[0]

                index = np.argmin(np.abs(bin_centers - n_used))
                temperature = temps[index]
                outputs = outputs / np.exp(temperature)

                outputs = softmax(outputs, axis=1)
                result[i, :] = outputs[0]

            avg_probs = np.mean(result, axis=0)
            best_idx = np.argmax(avg_probs)
            
            results_list.append({
                "sample": sample_name,
                "mixed class": class_names[best_idx],
                "mixed score": avg_probs[best_idx],
                "purity": "" 
            })
            
            print(f"Processed {filename} -> {class_names[best_idx]} ({avg_probs[best_idx]:.4f})")

        df_out = pd.DataFrame(results_list)
        df_out.to_csv(output_file, sep="\t", index=False)
        print(f"\nAll multi-results saved to {output_file}")

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
    parser_single.add_argument("-o", "--output_file", required=True, help="Path to save the output csv")
    parser_single.add_argument("-s", "--num_samplings", help="Number of subsamplings")
    parser_single.add_argument("-f", "--file_type", required=True, help="Type of file (e.g., bed, csv, bed_carlo)")

    # Parser for `sturgeon-v2 multi`
    parser_multi = subparsers.add_parser("multi", help="Predict on a directory of .bed or .csv files")
    parser_multi.add_argument("-i", "--input_dir", required=True, help="Input directory path")
    parser_multi.add_argument("-m", "--model", required=True, help="Path to the model zip")
    parser_multi.add_argument("-c", "--num_CpGs", help="Number of CpGs to use")
    parser_multi.add_argument("-o", "--output_file", required=True, help="Path to save the combined output tsv")
    parser_multi.add_argument("-s", "--num_samplings", help="Number of subsamplings")
    parser_multi.add_argument("-f", "--file_type", required=True, help="Type of files (e.g., bed, csv, bed_carlo)")

    args = parser.parse_args()
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
        )
    elif args.command == "multi":
        predict_directory(
            input_dir=args.input_dir,
            model_zip_path=args.model,
            n=args.num_CpGs,
            output_file=args.output_file,
            num_samples=num_samples,
            file_type=args.file_type,
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
