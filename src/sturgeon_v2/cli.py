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


# ---------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------
def softmax(x, axis=1):
    x = x - np.max(x, axis=axis, keepdims=True)
    exp_x = np.exp(x)
    return exp_x / np.sum(exp_x, axis=axis, keepdims=True)


# ---------------------------------------------------------------------
# Prediction
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

    # --------------------------------------------------------------
    # Extract model zip to temp directory
    # --------------------------------------------------------------
    tmpdir = tempfile.mkdtemp(prefix="sturgeon_onnx_")

    try:
        with zipfile.ZipFile(model_zip_path, "r") as zf:
            zf.extractall(tmpdir)

        model_root = os.path.join(tmpdir, model_zip_path.split('_')[-1].split('.')[0])
        model_path = os.path.join(model_root, "model.onnx")

        if not os.path.exists(model_path):
            raise FileNotFoundError(f"model.onnx not found at {model_path}")

        # --------------------------------------------------------------
        # Load classification config
        # --------------------------------------------------------------
        with open(os.path.join(model_root, "classification_system.yaml"), "r") as f:
            classification_file = yaml.safe_load(f)

        class_names = list(classification_file["encoder"]["type"].keys())
        output_size = len(class_names)

        # --------------------------------------------------------------
        # Load probe file
        # --------------------------------------------------------------
        probe_df = pl.read_csv(
            os.path.join(model_root, "probe.bed"),
            separator="\t",
        )

        input_size = len(probe_df)

        # --------------------------------------------------------------
        # Load input methylation data
        # --------------------------------------------------------------
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

        # --------------------------------------------------------------
        # Join probes with input
        # --------------------------------------------------------------
        nn_input = probe_df.join(
            bed_file,
            left_on="name",
            right_on="probe_id",
            how="left",
            validate="1:1",
        ).with_columns(pl.all().fill_null(0))

        nn_input = nn_input["methylation_call"].to_numpy().astype(np.float32)

        # --------------------------------------------------------------
        # Load ONNX model
        # --------------------------------------------------------------
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

        # --------------------------------------------------------------
        # Load calibration data
        # --------------------------------------------------------------
        bin_centers = np.load(os.path.join(model_root, "bin_centers.npy"))
        temps = np.load(os.path.join(model_root, "temps.npy"))

        # --------------------------------------------------------------
        # Run predictions
        # --------------------------------------------------------------
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

        # --------------------------------------------------------------
        # Write output
        # --------------------------------------------------------------
        df = pd.DataFrame(result, columns=class_names)
        df["probes"] = n_used
        df.to_csv(output_file, index=False)

    finally:
        # Cleanup temp directory AFTER inference
        shutil.rmtree(tmpdir, ignore_errors=True)


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input_file", required=True)
    parser.add_argument("-m", "--model", required=True)
    parser.add_argument("-c", "--num_CpGs")
    parser.add_argument("-o", "--output_file", required=True)
    parser.add_argument("-s", "--num_samplings")
    parser.add_argument("-f", "--file_type", required=True)

    args = parser.parse_args()
    num_samples = int(args.num_samplings) if args.num_samplings else 1

    print("---------------------------------------------------------------------")
    print("Running intraoperative methylation based classification (ONNX)")
    print("---------------------------------------------------------------------")

    predict(
        input_file=args.input_file,
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

