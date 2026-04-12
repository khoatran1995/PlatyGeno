import os
import torch
import numpy as np
import pandas as pd
from tqdm import tqdm
from platygeno import PlatyGenoEngine
from platygeno.evo_reader import read_evo_features

def analyze_activation_profile(input_path, engine, sample_size=1000):
    """
    Scans a small subset of reads to build an activation distribution profile.
    This helps find the 'Sweet Spot' threshold for novel gene discovery.
    """
    print("="*70)
    print(f"📊 Calibration Scan: {os.path.basename(input_path)}")
    print(f"🎯 Objective: Determine optimal discovery threshold")
    print("="*70)
    
    # 1. Run a 'Wide Net' scan with zero threshold
    # We use platygeno's internal reader but with no filters
    results = read_evo_features(input_path, engine, start=0, stop=sample_size)
    
    if results.empty:
        print("❌ Error: No features detected in sample subset.")
        return

    # 2. Analyze the 'Max Activation' per read
    max_activations = results.groupby('read_id')['activation'].max().values
    
    # 3. Calculate Percentiles
    stats = {
        "Mean": np.mean(max_activations),
        "Max": np.max(max_activations),
        "Median": np.median(max_activations),
        "Top 5%": np.percentile(max_activations, 95),
        "Top 1%": np.percentile(max_activations, 99),
        "Top 0.1%": np.percentile(max_activations, 99.9)
    }

    print("\n📈 Activation Distribution (Outlier Profile):")
    for label, val in stats.items():
        print(f"   {label:10}: {val:.4f}")

    print("\n💡 Recommendation for Discovery:")
    if stats["Top 1%"] > 8.0:
        print(f"   🔥 High Signal Sample: Set min_activation to 8.0 or 9.0")
    elif stats["Top 1%"] > 5.0:
        print(f"   🔎 Subtle Sample: Set min_activation to 6.5 or 7.0")
    else:
        print(f"   ⚠️ Low Signal Sample: Set min_activation to 4.0 or 5.0")
        print("      (Or download more reads to find rare extreme outliers)")
    print("="*70)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Calibrate PlatyGeno discovery thresholds.")
    parser.add_argument("--input", type=str, required=True, help="Path to sample FASTQ")
    parser.add_argument("--n", type=int, default=1000, help="Number of reads to scan for calibration")
    args = parser.parse_args()

    # Initialize Engine
    engine = PlatyGenoEngine()
    analyze_activation_profile(args.input, engine, sample_size=args.n)
