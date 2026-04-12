import os
import argparse
import pandas as pd
import platygeno

def run_discovery_showcase(input_path=None, top_pct=None, start=0, limit=None):
    print("="*70)
    print("PHASE 1: Iterative Genomic Discovery (5,000 Reads Mode)")
    print("="*70)
    
    # 1. Configuration
    if input_path is None:
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        input_path = os.path.join(base_dir, "data", "SRR23196177_subset.fastq")
    
    if not os.path.exists(input_path):
        print(f"Error: {input_path} not found.")
        return

    scan_end = start + limit if limit else None
    
    # 2. Sequential Discovery
    print(f"Scanning Reads: {start} to {scan_end if scan_end else 'End'}")
    
    # Dynamic filename for chunked runs
    output_label = f"{start}_{limit}" if limit else "full"
    output_csv = f"discovery_hits_{output_label}.csv"

    results = platygeno.discover_genes(
        input_path=input_path,
        scan_start=start,
        scan_end=scan_end, 
        min_activation=5.0, 
        rel_freq_max=0.001, # Adjusted to 0.1% (allowing up to 5 reads per 5k chunk)
        top_n=20, 
        top_pct=top_pct,
        output_path=output_csv
    )
    
    if results.empty:
        print(f"⚠️ Chunk {start}-{scan_end} yielded 0 novel features.")
        return

    print("\nCHUNK SCAN COMPLETE")
    print("="*70)
    print(f"Discovery Report: {output_csv}")
    print(f"Detected {len(results)//2} winners.")
    print(f"NEXT STEP: Run 'python validation/step2_local_blast.py'")
    print("="*70)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PlatyGeno: Chunked Discovery.")
    parser.add_argument("--input", type=str, help="Path to input FASTQ/FASTA file")
    parser.add_argument("--start", type=int, default=0, help="Read index to start scanning from (e.g. 0)")
    parser.add_argument("--limit", type=int, default=5000, help="Number of reads to scan (e.g. 5000)")
    parser.add_argument("--top-pct", type=float, help="Top % of outliers to target (e.g. 0.01)")
    
    args = parser.parse_args()
    
    run_discovery_showcase(
        input_path=args.input, 
        top_pct=args.top_pct,
        start=args.start,
        limit=args.limit
    )
