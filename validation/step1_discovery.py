import os
import argparse
import pandas as pd
import platygeno

# Known features to subtract by default (Digital Subtraction)
DEFAULT_SUPPRESSION = [212, 16509, 32214, 30446, 26886, 26280, 22717]

def run_discovery_showcase(input_path=None, top_pct=None, start=0, limit=None, exclude_ids=None):
    print("="*70)
    print("PHASE 1: Iterative Genomic Discovery (Scale-Aware + Suppression)")
    print("="*70)
    
    # 1. Configuration
    if input_path is None:
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        input_path = os.path.join(base_dir, "data", "SRR23196177_subset.fastq")
    
    if not os.path.exists(input_path):
        print(f"Error: {input_path} not found.")
        return

    # Combine default suppression with user-provided IDs
    exclusion_list = DEFAULT_SUPPRESSION.copy()
    if exclude_ids:
        exclusion_list.extend(exclude_ids)
    exclusion_list = list(set(exclusion_list)) # Unique

    scan_end = start + limit if limit else None
    
    # 2. Sequential Discovery
    # Using 'Deep-Search' settings (Threshold 5.0+, Rarity 0.01%)
    print(f"Scanning Reads: {start} to {scan_end if scan_end else 'End'}")
    print(f"Surgical Sensitivity: 5.0 | Rarity Filter: 0.01%")
    
    # Dynamic filename for chunked runs
    output_label = f"{start}_{limit}" if limit else "full"
    output_csv = f"discovery_hits_{output_label}.csv"

    results = platygeno.discover_genes(
        input_path=input_path,
        scan_start=start,
        scan_end=scan_end, 
        min_activation=5.0, 
        rel_freq_max=0.0001, # Stricter Rarity Filter for chunked runs
        top_n=20, 
        top_pct=top_pct,
        excluded_features=exclusion_list, # DIGITAL SUBTRACTION
        output_path=output_csv
    )
    
    if results.empty:
        print(f"⚠️ Chunk {start}-{scan_end} yielded 0 novel features.")
        return

    print("\nCHUNK SCAN COMPLETE")
    print("="*70)
    print(f"Discovery Report: {output_csv}")
    print(f"Detected {len(results)//2} winners (suppressing {len(exclusion_list)} defaults).")
    print(f"NEXT STEP: Run 'python validation/step2_local_blast.py' using {output_csv}")
    print("="*70)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PlatyGeno: Chunked Discovery with Digital Subtraction.")
    parser.add_argument("--input", type=str, help="Path to input FASTQ/FASTA file")
    parser.add_argument("--start", type=int, default=0, help="Read index to start scanning from (e.g. 0)")
    parser.add_argument("--limit", type=int, default=5000, help="Number of reads to scan (e.g. 5000)")
    parser.add_argument("--top-pct", type=float, help="Top % of outliers to target (e.g. 0.01)")
    parser.add_argument("--exclude", type=str, help="Comma-separated Feature IDs to subtract (e.g. 212,16509)")
    
    args = parser.parse_args()
    
    exclude_list = []
    if args.exclude:
        try:
            exclude_list = [int(x.strip()) for x in args.exclude.split(",")]
        except ValueError:
            print("Error: Exclusion list must be comma-separated integers.")

    run_discovery_showcase(
        input_path=args.input, 
        top_pct=args.top_pct,
        start=args.start,
        limit=args.limit,
        exclude_ids=exclude_list
    )
