#!/usr/bin/env python3

import argparse
import math
import sys
import os
from pathlib import Path
import pandas as pd
import torch
import warnings
from Bio import SeqIO

warnings.filterwarnings("ignore",
    message="Setting attributes on ParameterList is not supported.")


def load_model(seq2topt_dir, model_type="topt"):
    code_dir = os.path.join(seq2topt_dir, "code")
    if code_dir not in sys.path:
        sys.path.insert(0, code_dir)

    from model import MultiAttModel

    dim = 320
    window = 3
    n_head = 4
    n_RD = 4

    model = MultiAttModel(dim, window, n_head, n_RD)

    weights_dir = os.path.join(seq2topt_dir, "weights")
    if model_type == "topt":
        weight_file = os.path.join(weights_dir, "model_topt_window.3_r2.0.57.pth")
        scale_factor = 120.0
    elif model_type == "tm":
        weight_file = os.path.join(weights_dir, "model_tm_window.3_r2.0.76.pth")
        scale_factor = 100.0
    else:
        raise ValueError(f"Unknown model_type: {model_type}")

    if not os.path.isfile(weight_file):
        sys.exit(f"[ERROR] Weight file not found: {weight_file}")

    state_dict = torch.load(weight_file, map_location="cpu", weights_only=False)
    model.load_state_dict(state_dict)
    model.eval()

    return model, scale_factor


def get_esm_model():
    import esm
    model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
    model.eval()
    batch_converter = alphabet.get_batch_converter()
    return model, batch_converter


def predict_sequences(fasta_path, seq2topt_dir, model_type="topt"):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return pd.DataFrame(columns=["seq_id", "prediction"])

    print(f"[INFO] Loading ESM-2 model...", file=sys.stderr)
    esm_model, batch_converter = get_esm_model()

    print(f"[INFO] Loading Seq2{model_type.capitalize()} model...", file=sys.stderr)
    model, scale_factor = load_model(seq2topt_dir, model_type)

    seq_ids = [rec.id for rec in records]
    seq_strs = [str(rec.seq) for rec in records]

    predictions = []
    batch_size = 4

    for i in range(math.ceil(len(seq_ids) / batch_size)):
        batch_ids = seq_ids[i * batch_size: (i + 1) * batch_size]
        batch_seqs = seq_strs[i * batch_size: (i + 1) * batch_size]

        inputs = [(batch_ids[j], batch_seqs[j]) for j in range(len(batch_ids))]

        try:
            _, _, batch_tokens = batch_converter(inputs)

            with torch.no_grad():
                results_esm = esm_model(
                    batch_tokens, repr_layers=[6], return_contacts=False
                )
                emb = results_esm["representations"][6]
                emb = emb.transpose(1, 2)

                preds = model(emb)

            pred_values = preds.cpu().detach().numpy().reshape(-1).tolist()
            predictions.extend(pred_values)

        except Exception as e:
            print(f"[WARN] Batch {i} failed: {e}", file=sys.stderr)
            predictions.extend([float("nan")] * len(batch_ids))

    scaled = [v * scale_factor for v in predictions]

    return pd.DataFrame({
        "seq_id": seq_ids,
        "prediction": scaled,
    })


def main():
    p = argparse.ArgumentParser(
        description="Run Seq2Topt/Seq2Tm predictions on a FASTA file."
    )
    p.add_argument("--fasta", required=True,
                   help="Input multi-FASTA file")
    p.add_argument("--seq2topt_dir", required=True,
                   help="Path to cloned Seq2Topt repository")
    p.add_argument("--output_topt", default=None,
                   help="Output TSV for Topt predictions")
    p.add_argument("--output_tm", default=None,
                   help="Output TSV for Tm predictions")
    args = p.parse_args()

    if not args.output_topt and not args.output_tm:
        sys.exit("[ERROR] At least one of --output_topt or --output_tm required.")

    fasta = Path(args.fasta).resolve()
    seq2topt_dir = Path(args.seq2topt_dir).resolve()

    if not fasta.exists():
        sys.exit(f"[ERROR] FASTA file not found: {fasta}")

    if args.output_topt:
        print("[INFO] Predicting Topt...", file=sys.stderr)
        df_topt = predict_sequences(fasta, seq2topt_dir, model_type="topt")
        df_topt.columns = ["seq_id", "predicted_topt_C"]
        out = Path(args.output_topt)
        out.parent.mkdir(parents=True, exist_ok=True)
        df_topt.to_csv(out, sep="\t", index=False)
        print(f"[INFO] Wrote {len(df_topt)} Topt predictions to {out}",
              file=sys.stderr)

    if args.output_tm:
        print("[INFO] Predicting Tm...", file=sys.stderr)
        df_tm = predict_sequences(fasta, seq2topt_dir, model_type="tm")
        df_tm.columns = ["seq_id", "predicted_tm_C"]
        out = Path(args.output_tm)
        out.parent.mkdir(parents=True, exist_ok=True)
        df_tm.to_csv(out, sep="\t", index=False)
        print(f"[INFO] Wrote {len(df_tm)} Tm predictions to {out}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
