# Copyright 2026 Khoa Tu Tran
import argparse
import os
import sys
import platygeno

def main():
    parser = argparse.ArgumentParser(description="PlatyGeno: Unsupervised Gene Discovery via Evo 2 & SAE")
    
    # 1. Input/Output
    parser.add_argument("--input", "-i", type=str, required=True, help="Input FASTQ/FASTA file")
    parser.add_argument("--output", "-o", type=str, default="discovery_results.csv", help="Output results CSV")
    
    # 2. Scanning Range
    parser.add_argument("--start", type=int, default=0, help="First read index (default: 0)")
    parser.add_argument("--limit", "-l", type=int, default=5000, help="Number of reads to scan (default: 5000)")
    parser.add_argument("--batch-size", "-b", type=int, default=16, help="GPU batch size (default: 16)")
    
    # 3. Discovery Logic
    parser.add_argument("--threshold", "-t", type=float, default=5.0, help="Min activation (default: 5.0)")
    parser.add_argument("--rarity-only", action="store_true", help="Enable rarity filtering to target novel dark matter")
    parser.add_argument("--top-n", "-n", type=int, default=25, help="Number of features to target (default: 25)")
    parser.add_argument("--exclude", type=str, help="Comma-separated feature IDs to ignore (e.g. 212,16509)")
    
    # 4. Assembly 
    parser.add_argument("--window", "-w", type=int, default=60, help="Snippet window size (bp)")
    parser.add_argument("--overlap", type=int, default=20, help="Min assembly overlap (bp)")
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
        
    args = parser.parse_args()

    # Parse exclusions
    excluded = [int(x.strip()) for x in args.exclude.split(",")] if args.exclude else None
    
    # Logic: Panoramic (1.0) is now the default.
    rel_freq = 0.001 if args.rarity_only else 1.0

    # MASTER PIPELINE CALL
    results = platygeno.discover_genes(
        input_path=args.input,
        scan_start=args.start,
        scan_end=args.start + args.limit,
        top_n=args.top_n,
        window_size=args.window,
        min_overlap=args.overlap,
        min_activation=args.threshold,
        rel_freq_max=rel_freq,
        batch_size=args.batch_size,
        excluded_features=excluded,
        output_path=args.output
    )

    if not results.empty:
        print("\n🧬 DISCOVERY RESULTS SUMMARY")
        print("="*100)
        # Select key columns for the summary table
        cols = ['feature_id', 'feature_name', 'activation', 'rarity_pct', 'length']
        # Filter columns to only those that exist
        available_cols = [c for c in cols if c in results.columns]
        print(results[available_cols].to_string(index=False))
        print("="*100)
        print(f"✅ Success! Full report saved to {args.output}")
    else:
        print("\n⚠️ No genes discovered with the current parameters.")

if __name__ == "__main__":
    main()
