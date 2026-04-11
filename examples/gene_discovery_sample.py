# Copyright 2026 Khoa Tu Tran
import argparse
import os
import platygeno

def main():
    parser = argparse.ArgumentParser(description="PlatyGeno: End-to-End Gene Discovery")
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    # Path Arguments
    parser.add_argument("--input", type=str, default=os.path.join(base_dir, "data", "sample.fastq"),
                        help="Input FASTQ/FASTA file")
    parser.add_argument("--output", type=str, default=os.path.join(base_dir, "data", "gene_discovery_results.csv"),
                        help="Output results CSV")
    
    # Range & Filtering Arguments
    parser.add_argument("--start", type=int, default=0, help="First read index to scan")
    parser.add_argument("--end", type=int, default=4000, help="Last read index to scan")
    parser.add_argument("--threshold", type=float, default=5.0, help="Min activation to consider")
    
    # Extraction & Assembly Arguments
    parser.add_argument("--top_n", type=int, default=10, help="Number of rare features to target")
    parser.add_argument("--window", type=int, default=60, help="Snippet window size (bp)")
    parser.add_argument("--min_overlap", type=int, default=20, help="Min overlap for assembly (bp)")
    
    args = parser.parse_args()

    # The entire 3-phase workflow is now reduced to this single line:
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
        # Print a clean summary table for the user
        print("\n" + "="*60)
        print("🧬 DISCOVERY RESULTS SUMMARY")
        print("="*60)
        summary = results[['method', 'feature_id', 'activation', 'length']]
        print(summary.to_string(index=False))
        print("="*60)

if __name__ == "__main__":
    main()
