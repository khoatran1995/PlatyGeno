import os
import sys
import argparse
import subprocess

def run_suite(input_path, limit=20000, batch_size=32, threads=5, panoramic=True):
    print("="*80)
    print("      PLATYGENO DISCOVERY SUITE: 20K BENCHMARK PIPELINE")
    print("="*80)
    
    # 1.1 Dynamic Naming Prep
    base_name = os.path.splitext(os.path.basename(input_path))[0]
    csv_name = os.path.join("results", f"PLG_{base_name}_Significance.csv")

    # 1.2 Step 1: Significance Discovery
    discovery_cmd = [
        sys.executable, "validation/step1_discovery.py",
        "--input", input_path,
        "--limit", str(limit),
        "--batch-size", str(batch_size)
    ]
    
    print(f"\n🚀 STEP 1: Running significance scan on {limit} reads...")
    result = subprocess.run(discovery_cmd)
    
    if result.returncode != 0:
        print("❌ Step 1 failed. Aborting.")
        return
    
    if not os.path.exists(csv_name):
        print(f"❌ Error: {csv_name} was not generated.")
        return

    # 2. Step 2: Turbo-BLAST Validation
    # In benchmark mode, we validate BOTH (Consensus + Snippets) per your request.
    blast_cmd = [
        sys.executable, "validation/step2_blast.py",
        "--input", csv_name,
        "--threads", str(threads),
        "--all"
    ]
    
    print(f"\n📡 STEP 2: Validating novel candidates via Turbo-BLAST...")
    subprocess.run(blast_cmd)

    print("\n" + "="*80)
    print("✅ 20K BENCHMARK COMPLETE")
    print(f"Discovery Map:   {csv_name}")
    print(f"Validation Map:  results/PLG_{base_name}_Validation.csv")
    print(f"Dashboard:       reports/{base_name}_Dashboard.html")
    print("-" * 80)
    print("📈 NEXT STEPS:")
    print(f"AI Scan Complete. You can now review 'reports/{base_name}_Dashboard.html'.")
    print("To verify your 'Unknown' discoveries, we suggest manually extracting the")
    print("sequences and running them through structural tools like AlphaFold 2.")
    print("="*80)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PlatyGeno: Unified Discovery & Validation Pipeline.")
    parser.add_argument("--input", type=str, required=True, help="Path to raw sequence file")
    parser.add_argument("--limit", type=int, default=20000, help="Number of reads to scan")
    parser.add_argument("--batch-size", type=int, default=32, help="GPU batch size")
    parser.add_argument("--threads", type=int, default=5, help="Parallel BLAST threads")
    parser.add_argument("--panoramic", action="store_true", default=True, help="Run in Panoramic mode")
    
    args = parser.parse_args()
    run_suite(args.input, args.limit, args.batch_size, args.threads, args.panoramic)
