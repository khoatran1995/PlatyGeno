def print_runpod_instructions():
    print("="*70)
    print("🎯 STEP 4: AlphaFold 2 (Gold Standard) Execution on RunPod")
    print("="*70)
    
    # Path to the FASTA package created in Step 4
    base_dir = os.path.dirname(__file__)
    fasta_path = os.path.join(base_dir, "potential_novel_sequences.fasta")
    results_dir = os.path.join(base_dir, "alpha_results")
    
    print("🚀 TARGETING STRUCTURAL VALIDATION OF NOVEL GENES")
    print("-" * 50)
    print(f"📖 Input File:   {fasta_path}")
    print(f"📂 Results Dir:  {results_dir}")
    print("-" * 50)
    
    print("\n📦 REQUIRED SETUP (Run once on your RunPod terminal):")
    print("-" * 50)
    print("pip install torch")
    print("pip install \"colabfold[alphafold,openmm] @ git+https://github.com/sokrypton/ColabFold\"")
    print("-" * 50)
    
    print("\n🏛️ THE GOLD STANDARD RUN COMMAND:")
    print("-" * 50)
    print(f"colabfold_batch {fasta_path} {results_dir}")
    print("-" * 50)
    
    print("\n🏆 SCIENTIFIC REPORTING:")
    print("Once complete, your 3D PDB models and pLDDT confidence scores will be")
    print("safely stored in the 'validation/alpha_results/' directory.")
    print("="*70)

if __name__ == "__main__":
    print_runpod_instructions()
