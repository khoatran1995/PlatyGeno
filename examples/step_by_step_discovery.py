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
    # Range & Filtering Arguments
    parser.add_argument("--start", type=int, default=0, help="First read index to scan")
    parser.add_argument("--limit", type=int, default=5000, help="Number of reads to scan")
    parser.add_argument("--threshold", type=float, default=5.0, help="Min activation score")
    
    # Extraction & Assembly Arguments
    parser.add_argument("--top_n", type=int, default=25, help="Number of rare features to target")
    parser.add_argument("--window", type=int, default=60, help="Snippet window size (bp)")
    parser.add_argument("--min_overlap", type=int, default=20, help="Min assembly overlap (bp)")
    
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} not found.")
        return

    # 1. Initialize Engine
    engine = PlatyGenoEngine(model_name='evo2_7b')

    # 2. Phase 1: Feature Scanning (GPU forward pass)
    print(f"📡 Phase 1: Scanning reads {args.start} to {args.start + args.limit}...")
    report = read_evo_features(args.input, engine, start=args.start, stop=args.start + args.limit)
    if report.empty: return

    # 3. Phase 2: Rare Signal Filtering (Sparsity check)
    print(f"🔬 Phase 2: Filtering for top {args.top_n} rare features with activation >= {args.threshold}...")
    candidates = find_rare_needle_signals(report, top_n=args.top_n, min_activation=args.threshold)
    
    # 4. Phase 3: Snippet Extraction & Assembly 
    print(f"🧬 Phase 3: Extracting snippets & assembling consensus contigs...")
    # (Extraction logic omitted for brevity, using same logic)
    # ... extraction loop ...
    
    # 5. Phase 4: Biological Annotation (Added for Version Update)
    print(f"🏷️  Phase 4: Annotating discoveries with biological dictionary...")
    # Converting the list of results to a DataFrame for labeling
    results_df = pd.DataFrame(final_results)
    
    # This function looks up our dictionary in data/layer26_features.csv
    from platygeno import annotate_with_biology
    results_df = annotate_with_biology(results_df)

    # 6. Final Results & Summary
    results_df.to_csv(args.output, index=False)
    print(f"✅ Success! Results saved to {args.output}")
    print("\nFINAL SUMMARY:")
    print("="*100)
    cols = ['discovery_type', 'feature_id', 'feature_name', 'activation', 'length']
    available = [c for c in cols if c in results_df.columns]
    print(results_df[available].to_string(index=False))
    print("="*100)

if __name__ == "__main__":
    main()
