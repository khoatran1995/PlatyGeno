import os
import pandas as pd

def translate_dna(dna):
    """Translates DNA to protein across all 6 frames."""
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    def get_longest_orf(seq):
        tokens = [seq[i:i+3] for i in range(0, (len(seq)//3)*3, 3)]
        protein = "".join([table.get(t, 'X') for t in tokens])
        orfs = protein.split('_')
        return max(orfs, key=len) if orfs else ""

    rev_comp = "".join([{'A':'T','C':'G','G':'C','T':'A'}.get(b,'N') for b in reversed(dna)])
    frames = [dna, dna[1:], dna[2:], rev_comp, rev_comp[1:], rev_comp[2:]]
    return max([get_longest_orf(f) for f in frames], key=len)

def prepare_alphafold_package():
    print("="*70)
    print("👑 PHASE 4: AlphaFold 2 Preparation (Gold Standard)")
    print("="*70)
    
    # Input CSV from discovery/validation phase
    input_csv = "ibd_clinical_report.csv"
    if not os.path.exists(input_csv):
        if os.path.exists("validation/ibd_clinical_report.csv"):
            input_csv = "validation/ibd_clinical_report.csv"
        else:
            print(f"❌ Error: {input_csv} not found. Please run Step 2/3 first.")
            return

    df = pd.read_csv(input_csv, encoding='utf-8-sig')
    print(f"✅ Loaded {len(df)} candidates for AlphaFold 2 validation.")

    # Create FASTA file for ColabFold
    fasta_path = "validation/alphafold_input.fasta"
    os.makedirs("validation/alpha_results", exist_ok=True)
    
    with open(fasta_path, "w") as f:
        for i, row in df.iterrows():
            dna = row['sequence']
            fid = row['feature_id']
            protein = translate_dna(dna)
            
            # Skip very short sequences (< 5 aa) for AlphaFold
            if len(protein) < 5:
                continue
                
            f.write(f">feature_{fid}\n{protein}\n")
    
    print(f"✅ Standardized FASTA package created at: {fasta_path}")
    print("\n🚀 TO RUN ALPHAFOLD 2 (COLABFOLD) ON YOUR POD:")
    print("-" * 50)
    print("1. Install ColabFold (One-time):")
    print("   pip install 'colabfold[alphafold,openmm] @ git+https://github.com/sokrypton/ColabFold'")
    print("\n2. Launch AlphaFold Batch Run:")
    print(f"   colabfold_batch {fasta_path} validation/alpha_results/")
    print("-" * 50)
    print("\nAlphaFold will now begin high-precision folding. This may take 1-2 hours.")

if __name__ == "__main__":
    prepare_alphafold_package()
