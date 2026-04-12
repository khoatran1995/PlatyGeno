import os
import argparse
import pandas as pd
import platygeno

def run_discovery_showcase(input_path=None, top_pct=None):
    print("="*70)
    print("PHASE 1: Rare Genomic Feature Discovery (Clinical Benchmark)")
    print("="*70)
    
    # 1. Configuration
    if input_path is None:
        # Default to the shared benchmark dataset
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        input_path = os.path.join(base_dir, "data", "SRR23196177_subset.fastq")
    
    if not os.path.exists(input_path):
        print(f"Error: {input_path} not found.")
        print("Please ensure you have downloaded the sample data first.")
        return

    # 2. Scanning (GPU Phase)
    # We use Surgical Sensitivity (5.0 threshold) and Scale-Aware Rarity (0.1%)
    print(f"Analyzing Metagenome with Surgical Sensitivity Mode (Threshold: 5.0)...")
    results = platygeno.discover_genes(
        input_path=input_path,
        scan_end=None, 
        min_activation=5.0, 
        rel_freq_max=0.001, # Scale-Aware rarity limit (0.1%)
        top_n=200,          # Fallback if top_pct is None
        top_pct=top_pct     # Dynamic target (e.g., top 1% of outliers)
    )
    
    if results.empty:
        print("⚠️ No rare features met the activation threshold (>=5.0).")
        print("Warning: No unique genomic features detected.")
        return

    # 3. Save Results
    output_csv = "discovery_hits.csv"
    results.to_csv(output_csv, index=False)
    
    print("\nSTEP 1 COMPLETE")
    print("="*70)
    print(f"Discovery Report: {output_csv}")
    print(f"Detected {len(results)//2} winning features ({len(results)} rows).")
    print("NEXT STEP: Run 'python validation/step2_local_blast.py'")
    print("="*70)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PlatyGeno Step 1: Discover rare genomic features.")
    parser.add_argument("--input", type=str, help="Path to input FASTQ/FASTA file")
    parser.add_argument("--rarity-pct", type=float, default=0.001, help="Relative rarity threshold (default 0.1%)")
    parser.add_argument("--top-pct", type=float, help="Top percentage of outliers to target (e.g. 0.01 for 1%)")
    args = parser.parse_args()
    
    # Execute the discovery showcase
    run_discovery_showcase(input_path=args.input, top_pct=args.top_pct)
