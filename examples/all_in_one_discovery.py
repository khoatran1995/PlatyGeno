# Copyright 2026 Khoa Tu Tran
# gene_discovery_pipeline.py
#
# A high-level, end-to-end example using the PlatyGeno Discovery Pipeline.
# This script uses the simple `platygeno.discover_genes()` API to 
# run the entire discovery process in a single function call.

import argparse
import os
import platygeno

def main():
    parser = argparse.ArgumentParser(description="PlatyGeno One-Line Discovery Pipeline")
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    # Simple Arguments
    parser.add_argument("--input", type=str, default=os.path.join(base_dir, "data", "sample.fastq"),
                        help="Input FASTQ/FASTA file")
    parser.add_argument("--output", type=str, default=os.path.join(base_dir, "data", "gene_discovery_results.csv"),
                        help="Output results CSV")
    
    # Range & Filtering Arguments
    parser.add_argument("--start", type=int, default=0, help="First read index to scan")
    parser.add_argument("--limit", type=int, default=5000, help="Number of reads to scan")
    parser.add_argument("--threshold", type=float, default=5.0, help="Min activation score")
    
    # Extraction & Assembly Arguments
    parser.add_argument("--top_n", type=int, default=25, help="Number of rare features to target")
    parser.add_argument("--window", type=int, default=60, help="Snippet window size (bp)")
    parser.add_argument("--min_overlap", type=int, default=20, help="Min assembly overlap (bp)")
    
    args = parser.parse_args()

    # MASTER PIPELINE CALL (Unified API)
    results = platygeno.discover_genes(
        input_path=args.input,
        scan_start=args.start,
        scan_end=args.start + args.limit,
        top_n=args.top_n,
        window_size=args.window,
        min_overlap=args.min_overlap,
        min_activation=args.threshold,
        output_path=args.output
    )

    if not results.empty:
        print("\n🧬 DISCOVERY RESULTS SUMMARY")
        print("="*100)
        # We now include the 'feature_name' (e.g. Found: Coding Region vs Unknown)
        summary = results[['discovery_type', 'feature_id', 'feature_name', 'activation', 'length']]
        print(summary.to_string(index=False))
        print("="*100)
        print(f"✅ Success! Full report saved to {args.output}")

if __name__ == "__main__":
    main()
