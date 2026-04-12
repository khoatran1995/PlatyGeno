import os
import argparse
import pandas as pd
import platygeno

def run_significance_scan(input_path=None, top_pct=None, start=0, limit=None, ignore_rarity=True, batch_size=16):
    print("="*70)
    print("PHASE 1: Reference-Free Significance Mapping (Bio-Beacon)")
    if ignore_rarity:
        print(f"MODE: Landmarks Mode (Zero-Reference | Batch Size: {batch_size})")
    else:
        print(f"MODE: Rarity Mode (Deep-Mining | Batch Size: {batch_size})")
    print("="*70)
    
    # 1. Configuration
    if input_path is None:
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        input_path = os.path.join(base_dir, "data", "SRR23196177_subset.fastq")
    
    if not os.path.exists(input_path):
        print(f"Error: {input_path} not found.")
        return

    scan_end = start + limit if limit else None
    
    # Lead with significance as the primary scientific driver
    min_activation = 3.0
    rel_freq_max = 1.0 if ignore_rarity else 0.001
    top_n = 50 if ignore_rarity else 20

    # 2. Bio-Beacon Discovery
    print(f"📡 Scanning for biological landmarks in reads {start} to {scan_end if scan_end else 'End'}...")
    
    # Scientific labeling for output
    mode_label = "panoramic" if ignore_rarity else "novelty"
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
        batch_size=batch_size,
        output_path=output_csv
    )
    
    if results.empty:
        print(f"⚠️ No high-confidence landmarks detected in this region.")
        return

    print(f"\nSIGNIFICANCE MAPPING COMPLETE")
    print("="*70)
    print(f"Landmark Report: {output_csv}")
    print(f"Mapped {len(results)//2} significant biological features.")
    print(f"NEXT STEP: Run 'python validation/step2_local_blast.py --input {output_csv}'")
    print("="*70)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PlatyGeno: Unsupervised Genomic Significance Scanner.")
    parser.add_argument("--input", type=str, help="Path to raw sequence file (FASTQ/FASTA)")
    parser.add_argument("--start", type=int, default=0, help="First read index to process")
    parser.add_argument("--limit", type=int, default=5000, help="Number of reads to scan")
    parser.add_argument("--top-pct", type=float, help="Select the top X% of significant hits")
    parser.add_argument("--batch-size", type=int, default=16, help="Number of sequences per GPU batch (default: 16)")
    parser.add_argument("--rarity-only", action="store_true", help="Enable rarity filtering to target novel dark matter")
    
    args = parser.parse_args()
    
    # Default is now Significance Scanning (ignore_rarity=True)
    run_significance_scan(
        input_path=args.input, 
        top_pct=args.top_pct,
        start=args.start,
        limit=args.limit,
        ignore_rarity=not args.rarity_only,
        batch_size=args.batch_size
    )
