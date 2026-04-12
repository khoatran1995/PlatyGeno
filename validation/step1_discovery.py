import os
import pandas as pd
import platygeno
import argparse

def run_discovery_showcase(input_path=None):
    print("="*70)
    print("Starting PlatyGeno GPU Discovery Phase (High-Speed)")
    
    # 1. Handle Input Path
    if input_path is None:
        input_path = "data/HSMA33OT_R1.fastq"
        
    if not os.path.exists(input_path):
        if os.path.exists("sample.fastq"):
            input_path = "sample.fastq"
        elif os.path.exists("data/raw/sample.fastq"):
            input_path = "data/raw/sample.fastq"
        else:
            print(f"Error: {input_path} not found.")
            return

    print(f"Sample: {os.path.basename(input_path)}")
    print("="*70)

    # 2. Scanning (GPU Phase)
    print(f"Analyzing Metagenome with Ultra-Sensitive Mode (Threshold: 3.0)...")
    results = platygeno.discover_genes(
        input_path=input_path,
        scan_end=None, 
        min_activation=3.0, 
        top_n=200 
    )
    
    if results.empty:
        print("Warning: No unique genomic features detected.")
        return

    # 3. Saving Hits
    print(f"Isolated {len(results)} high-confidence features.")
    
    # Save results in the validation folder
    output_path = os.path.join(os.path.dirname(__file__), "discovery_hits.csv")
    results.to_csv(output_path, index=False)
    
    print("\nSTEP 1 COMPLETE")
    print("="*70)
    print(f"Results saved to: {output_path}")
    print("NEXT STEP: Run 'python validation/step2_local_blast.py'")
    print("="*70)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PlatyGeno Discovery Phase")
    parser.add_argument("--input", type=str, help="Path to input FASTQ/FASTA file")
    args = parser.parse_args()
    
    run_discovery_showcase(input_path=args.input)
