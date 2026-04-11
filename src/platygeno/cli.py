# Copyright 2026 Khoa Tu Tran
import argparse
import os
import sys
import platygeno

def main():
    parser = argparse.ArgumentParser(description="PlatyGeno: Unsupervised Gene Discovery via Evo 2 & SAE")
    
    # Simple Arguments
    parser.add_argument("--input", "-i", type=str, required=True,
                        help="Input FASTQ/FASTA file")
    parser.add_argument("--output", "-o", type=str, default="gene_discovery_results.csv",
                        help="Output results CSV (default: gene_discovery_results.csv)")
    
    # Range & Filtering Arguments
    parser.add_argument("--start", type=int, default=0, help="First read index to scan (default: 0)")
    parser.add_argument("--end", type=int, default=4000, help="Last read index to scan (default: 4000)")
    parser.add_argument("--threshold", "-t", type=float, default=5.0, help="Min activation score (default: 5.0)")
    
    # Extraction & Assembly Arguments
    parser.add_argument("--top-n", type=int, default=10, help="Number of rare features to target (default: 10)")
    parser.add_argument("--window", "-w", type=int, default=60, help="Snippet window size (bp) (default: 60)")
    parser.add_argument("--min-overlap", type=int, default=20, help="Min assembly overlap (bp) (default: 20)")
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
        
    args = parser.parse_args()

    # MASTER PIPELINE CALL
    results = platygeno.discover_genes(
        input_path=args.input,
        scan_start=args.start,
        scan_end=args.end,
        top_n=args.top_n,
        window_size=args.window,
        min_overlap=args.min_overlap,
        min_activation=args.threshold,
        output_path=args.output
    )

    if not results.empty:
        print("\n🧬 DISCOVERY RESULTS SUMMARY")
        print("="*60)
        summary = results[['method', 'feature_id', 'activation', 'length']]
        print(summary.to_string(index=False))
        print("="*60)
        print(f"✅ Success! Full report saved to {args.output}")
    else:
        print("\n⚠️ No genes discovered with the current parameters.")

if __name__ == "__main__":
    main()
