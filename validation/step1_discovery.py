import os
import pandas as pd
import platygeno

def run_discovery_showcase():
    print("="*70)
    print("🧬 Starting PlatyGeno GPU Discovery Phase (High-Speed)")
    print("📊 Sample: HSMA33OT_R1 (Clinical IBD Metagenome)")
    print("="*70)
    
    input_path = "data/sample.fastq"
    
    # Check data/sample.fastq first
    if not os.path.exists(input_path):
        # Fallback to current dir
        if os.path.exists("sample.fastq"):
            input_path = "sample.fastq"
        elif os.path.exists("data/raw/sample.fastq"):
            input_path = "data/raw/sample.fastq"
        else:
            print("❌ Error: sample.fastq not found.")
            return

    # 1. Scanning (GPU Phase)
    print(f"🔍 Analyzing Clinical Metagenome with Ultra-Sensitive Mode (Threshold: 3.0)...")
    results = platygeno.discover_genes(
        input_path=input_path,
        scan_end=None, 
        min_activation=3.0, # Lowered for high-sensitivity novel gene detection
        top_n=200 # Increased to capture more low-activation candidates
    )
    
    if results.empty:
        print("⚠️ No unique genomic features detected.")
        return

    # 2. Saving Hits (Standardized Validation)
    print(f"✅ Isolated {len(results)} high-confidence clinical features.")
    
    # Save the raw DNA hits directly in the validation folder
    output_path = os.path.join(os.path.dirname(__file__), "discovery_hits.csv")
    results.to_csv(output_path, index=False)
    
    print("\n🏆 STEP 1 COMPLETE")
    print("="*70)
    print(f"Results saved to: {output_path}")
    print("🚀 NEXT STEP: Run 'python validation/step2_local_blast.py'")
    print("="*70)

if __name__ == "__main__":
    run_discovery_showcase()
