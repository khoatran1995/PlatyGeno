# Copyright 2026 Khoa Tu Tran
# gene_discovery_sample.py
#
# A detailed, step-by-step example of the PlatyGeno discovery workflow.
# Use this script if you want to see exactly how the individual phases 
# (Scan -> Filter -> Extract -> Assemble) work under the hood.

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
    extract_precise_gene_code,
    assemble_feature_consensus
)

def main():
    parser = argparse.ArgumentParser(description="PlatyGeno Detailed Discovery Workflow")
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    parser.add_argument("--input", type=str, default=os.path.join(base_dir, "data", "sample.fastq"))
    parser.add_argument("--output", type=str, default=os.path.join(base_dir, "data", "gene_discovery_results.csv"))
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--end", type=int, default=4000)
    parser.add_argument("--window", type=int, default=60)
    parser.add_argument("--min_overlap", type=int, default=20)
    parser.add_argument("--threshold", type=float, default=5.0)
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} not found.")
        return

    # 1. Initialize Engine
    engine = PlatyGenoEngine(model_name='evo2_7b')

    # 2. Phase 1: Feature Scanning
    print(f"📡 Phase 1: Scanning reads {args.start} to {args.end}...")
    report = read_evo_features(args.input, engine, start=args.start, stop=args.end)
    if report.empty: return

    # 3. Phase 2: Rare Signal Filtering
    # We find features that appear rarely but with very high activation scores.
    print(f"🔬 Phase 2: Filtering for rare features with activation >= {args.threshold}...")
    candidates = find_rare_needle_signals(report, top_n=10, min_activation=args.threshold)
    
    # 4. Phase 3: Extraction (Dual Mode)
    print(f"🧬 Phase 3: Extracting snippets & assembling contigs...")
    ext = os.path.splitext(args.input)[1].lower()
    fmt = "fastq" if ext in [".fastq", ".fq"] else "fasta"
    
    # Pre-collect sequences needed for winners
    all_winning_read_ids = set()
    for fid in candidates['feature_id']:
        f_reads = report[report['feature_id'] == fid].sort_values('activation', ascending=False)
        all_winning_read_ids.update(f_reads.head(10)['read_id'].tolist())

    seq_map = {rec.id.replace("/","_").replace("|","_").replace(":","_"): str(rec.seq) 
               for rec in SeqIO.parse(args.input, fmt) if rec.id.replace("/","_").replace("|","_").replace(":","_") in all_winning_read_ids}

    final_results = []

    for fid in candidates['feature_id']:
        f_reads = report[report['feature_id'] == fid].sort_values('activation', ascending=False)
        best_row = f_reads.iloc[0]
        
        # --- METHOD A: Single Best Snippet (The 60bp "Anchor") ---
        read_id = best_row['read_id']
        dna = seq_map.get(read_id)
        if dna:
            snippet, act, pos = extract_precise_gene_code(engine, dna, fid, window_size=args.window)
            final_results.append({
                "method": "Best Snippet",
                "feature_id": fid,
                "read_id": read_id,
                "activation": round(act, 4),
                "length": len(snippet),
                "sequence": snippet
            })

        # --- METHOD B: Feature-Centric Assembly (Reconstructing context) ---
        top_10_reads = f_reads.head(10)['read_id'].tolist()
        sequences_to_assemble = [seq_map[rid] for rid in top_10_reads if rid in seq_map]
        
        if len(sequences_to_assemble) > 1:
            contig = assemble_feature_consensus(sequences_to_assemble, min_overlap=args.min_overlap)
            final_results.append({
                "method": "Assembled Contig",
                "feature_id": fid,
                "read_id": "Consensus (N=" + str(len(sequences_to_assemble)) + ")",
                "activation": round(best_row['activation'], 4),
                "length": len(contig),
                "sequence": contig
            })

    # Save Results
    results_df = pd.DataFrame(final_results)
    results_df.to_csv(args.output, index=False)
    print(f"✅ Success! Results saved to {args.output}")
    print(results_df[['method', 'feature_id', 'activation', 'length']].to_string(index=False))

if __name__ == "__main__":
    main()
