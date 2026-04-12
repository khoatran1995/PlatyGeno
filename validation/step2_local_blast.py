import os
import pandas as pd
import platygeno

def run_benchmarking_validation():
    """
    Standard Ph.D. Workflow: Step 2 - Novelty Validation
    This script is now a simple wrapper around the platygeno engine's 
    built-in validator.
    """
    input_csv = "discovery_hits.csv"
    output_csv = "blast_results.csv"
    novel_csv = "potential_novel_sequences.csv"

    if not os.path.exists(input_csv):
        print(f"❌ Error: {input_csv} not found. Please run Step 1 first.")
        return

    # Load discovery hits
    df = pd.read_csv(input_csv)
    
    # Execute native library validation
    # (The validator handles checkpointing and incremental saving internally)
    final_df = platygeno.validate_novelty(
        df, 
        output_path=output_csv, 
        novel_path=novel_csv
    )

    print("\n" + "="*70)
    print("STEP 2 COMPLETE: Novelty Validation")
    print("="*70)
    print(f"Audit Trail: {output_csv}")
    print(f"Candidates:  {novel_csv}")
    print("NEXT STEP: Run 'python validation/step3_fasta_prep.py'")
    print("="*70)

if __name__ == "__main__":
    run_benchmarking_validation()
