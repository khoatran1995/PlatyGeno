import os
import argparse
import pandas as pd
import platygeno

def run_discovery_showcase(input_path=None, top_pct=None, start=0, limit=None, is_panoramic=False):
    print("="*70)
    print("PHASE 1: Iterative Genomic Discovery (Scale-Aware)")
    if is_panoramic:
        print("MODE: Panoramic Mode (Disabling Rarity Filter - Showing EVERYTHING >3.0)")
    print("="*70)
    
    # 1. Configuration
    if input_path is None:
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        input_path = os.path.join(base_dir, "data", "SRR23196177_subset.fastq")
    
    if not os.path.exists(input_path):
        print(f"Error: {input_path} not found.")
        return

    scan_end = start + limit if limit else None
    
    # Define settings based on mode
    if is_panoramic:
        min_activation = 3.0
        rel_freq_max = 1.0  # 100% rarity allowed (Panasonic View)
        top_n = 100         # Allow a large winner circle
    else:
        min_activation = 5.0
        rel_freq_max = 0.001 # 0.1% Rarity Filter
        top_n = 20

    # 2. Sequential Discovery
    print(f"Scanning Reads: {start} to {scan_end if scan_end else 'End'}")
    
    # Dynamic filename for chunked runs
    mode_label = "panoramic" if is_panoramic else "surgical"
    output_label = f"{start}_{limit}_{mode_label}" if limit else mode_label
    output_csv = f"discovery_hits_{output_label}.csv"

    results = platygeno.discover_genes(
        input_path=input_path,
        scan_start=start,
        scan_end=scan_end, 
        min_activation=min_activation, 
        rel_freq_max=rel_freq_max, 
        top_n=top_n, 
        top_pct=top_pct,
        output_path=output_csv
    )
    
    if results.empty:
        print(f"⚠️ Chunk {start}-{scan_end} yielded 0 features.")
        return

    print(f"\n{mode_label.upper()} SCAN COMPLETE")
    print("="*70)
    print(f"Discovery Report: {output_csv}")
    print(f"Detected {len(results)//2} winning features.")
    print(f"NEXT STEP: Run 'python validation/step2_local_blast.py --input {output_csv}'")
    print("="*70)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PlatyGeno: Iterative Discovery.")
    parser.add_argument("--input", type=str, help="Path to input FASTQ/FASTA file")
    parser.add_argument("--start", type=int, default=0, help="Read index to start scanning (e.g. 0)")
    parser.add_argument("--limit", type=int, default=5000, help="Number of reads to scan (e.g. 5000)")
    parser.add_argument("--top-pct", type=float, help="Top % of outliers to target (e.g. 0.01)")
    parser.add_argument("--panoramic", action="store_true", help="Disable rarity filter and show all high-energy hits")
    
    args = parser.parse_args()
    
    run_discovery_showcase(
        input_path=args.input, 
        top_pct=args.top_pct,
        start=args.start,
        limit=args.limit,
        is_panoramic=args.panoramic
    )
