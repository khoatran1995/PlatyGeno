import os
import pandas as pd

def translate_dna(dna):
    """Simple DNA to Protein translation."""
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein = ""
    for i in range(0, len(dna)-(len(dna)%3), 3):
        codon = dna[i:i+3].upper()
        protein += table.get(codon, 'X')
    return protein

def prepare_fasta():
    print("="*70)
    print("🧬 STEP 3: Potential Novel Sequence FASTA Preparation")
    print("="*70)
    
    # Path to the novel hits generated in Step 2
    base_dir = os.path.dirname(__file__)
    input_csv = os.path.join(base_dir, "potential_novel_sequences.csv")
    output_fasta = os.path.join(base_dir, "potential_novel_sequences.fasta")
    
    if not os.path.exists(input_csv):
        print(f"❌ Error: {input_csv} not found. Please run Step 2 first.")
        return

    df = pd.read_csv(input_csv, encoding='utf-8-sig')
    print(f"✅ Loaded {len(df)} potential novel features.")

    # Create FASTA package
    with open(output_fasta, "w") as f:
        for i, row in df.iterrows():
            dna = row['sequence']
            fid = row['feature_id']
            protein = translate_dna(dna)
            
            # Skip short sequences or stop-codon interrupted sequences for folding
            if len(protein) < 5 or "_" in protein:
                continue
                
            f.write(f">feature_{fid}\n{protein}\n")
    
    print(f"✅ Standardized FASTA package created: {output_fasta}")
    print("\n🏆 STEP 3 COMPLETE")
    print("="*70)
    print("🚀 NEXT STEP: Run 'python validation/step4_alphafold_run.py'")
    print("="*70)

if __name__ == "__main__":
    prepare_fasta()
