# Copyright 2026 Khoa Tu Tran
import argparse
import os
import torch
import _codecs
import pandas as pd
from Bio import SeqIO

# Torch serialization safety bypass
torch.serialization.add_safe_globals([_codecs.encode])

from platygeno import (
    PlatyGenoEngine, 
    read_evo_features, 
    find_rare_needle_signals, 
    get_best_reads_for_features, 
    extract_precise_gene_code
)

def main():
    parser = argparse.ArgumentParser(description="PlatyGeno Gene Discovery Example")
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    parser.add_argument("--input", type=str, default=os.path.join(base_dir, "data", "sample.fastq"))
    parser.add_argument("--output", type=str, default=os.path.join(base_dir, "data", "gene_snippets.csv"))
    parser.add_argument("--limit", type=int, default=4000)
    parser.add_argument("--window", type=int, default=60)
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} not found.")
        return

    # Initialize Engine
    engine = PlatyGenoEngine(model_name='evo2_7b')

    # Phase 1: Feature Scanning
    print(f"📡 Scanning {args.limit} reads...")
    report = read_evo_features(args.input, engine, limit=args.limit)
    if report.empty: return

    # Phase 2: Filtering
    candidates = find_rare_needle_signals(report, top_n=10)
    winners = get_best_reads_for_features(report, candidates['feature_id'].tolist())

    # Phase 3: Extraction
    print(f"🧬 Extracting snippets (Window: {args.window}bp)...")
    winning_ids = set(winners['read_id'].tolist())
    ext = os.path.splitext(args.input)[1].lower()
    fmt = "fastq" if ext in [".fastq", ".fq"] else "fasta"
    
    # Map read_id back to sequence
    seq_map = {rec.id.replace("/","_").replace("|","_").replace(":","_"): str(rec.seq) 
               for rec in SeqIO.parse(args.input, fmt) if rec.id.replace("/","_").replace("|","_").replace(":","_") in winning_ids}

    final_results = []
    for _, row in winners.iterrows():
        feat_id, read_id = int(row['feature_id']), row['read_id']
        dna = seq_map.get(read_id)
        if not dna: continue

        snippet, act, pos = extract_precise_gene_code(engine, dna, feat_id, window_size=args.window)
        final_results.append({
            "feature_id": feat_id,
            "read_id": read_id,
            "activation": round(act, 4),
            "peak_pos": pos,
            "snippet": snippet
        })

    # Save Results
    results_df = pd.DataFrame(final_results)
    results_df.to_csv(args.output, index=False)
    print(f"✅ Success! Results saved to {args.output}")
    print(results_df.to_string(index=False))

if __name__ == "__main__":
    main()
