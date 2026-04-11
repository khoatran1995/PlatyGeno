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
    print(f"🔍 Analyzing Clinical Metagenome ({os.path.getsize(input_path)/1024:.1f} KB)...")
    results = platygeno.discover_genes(
        input_path=input_path,
        scan_end=None, 
        min_activation=5.0,
        top_n=15 # Increased for better discovery
    )
    
    if results.empty:
        print("⚠️ No unique genomic features detected.")
        return

    # 2. Saving Hits (Fast Exit)
    print(f"✅ Isolated {len(results)} high-confidence clinical features.")
    
    # Save the raw DNA hits for local validation
    results.to_csv("discovery_hits.csv", index=False)
    
    print("\n🏆 DISCOVERY PHASE COMPLETE")
    print("="*70)
    print(f"Discovery Results saved to: discovery_hits.csv")
    print("🚀 ACTION: Download 'discovery_hits.csv' and shut down your RunPod to save costs!")
    print("Then run 'python examples/validate_local.py' on your personal computer.")
    print("="*70)

if __name__ == "__main__":
    run_discovery_showcase()
