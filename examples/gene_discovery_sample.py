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
    extract_precise_gene_code,
    assemble_feature_consensus
)

def main():
    parser = argparse.ArgumentParser(description="PlatyGeno Gene Discovery Example")
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    parser.add_argument("--input", type=str, default=os.path.join(base_dir, "data", "sample.fastq"))
    parser.add_argument("--output", type=str, default=os.path.join(base_dir, "data", "gene_discovery_results.csv"))
    parser.add_argument("--limit", type=int, default=4000)
    parser.add_argument("--window", type=int, default=60)
    parser.add_argument("--min_overlap", type=int, default=20)
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

    # Phase 2: Filtering for Rare/High-Activation Candidates
    candidates = find_rare_needle_signals(report, top_n=10)
    
    print(f"🧬 Extracting snippets & assembling contigs (Window: {args.window}bp)...")
    ext = os.path.splitext(args.input)[1].lower()
    fmt = "fastq" if ext in [".fastq", ".fq"] else "fasta"
    
    # Pre-collect sequences needed for winners (optimization)
    all_winning_read_ids = set()
    for fid in candidates['feature_id']:
        f_reads = report[report['feature_id'] == fid].sort_values('activation', ascending=False)
        all_winning_read_ids.update(f_reads.head(10)['read_id'].tolist())

    seq_map = {rec.id.replace("/","_").replace("|","_").replace(":","_"): str(rec.seq) 
               for rec in SeqIO.parse(args.input, fmt) if rec.id.replace("/","_").replace("|","_").replace(":","_") in all_winning_read_ids}

    final_results = []

    for fid in candidates['feature_id']:
        # Get reads for this specific feature
        f_reads = report[report['feature_id'] == fid].sort_values('activation', ascending=False)
        
        # --- METHOD A: Single Best Snippet (Old Method) ---
        best_row = f_reads.iloc[0]
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

        # --- METHOD B: Feature-Centric Assembly (New Method) ---
        # Collect top 10 sequences for this feature
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
    
    # Print a clean summary table
    summary = results_df[['method', 'feature_id', 'activation', 'length']]
    print(summary.to_string(index=False))

if __name__ == "__main__":
    main()
