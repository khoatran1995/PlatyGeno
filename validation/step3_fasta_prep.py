import os
import argparse
import pandas as pd

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(dna))

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

def find_longest_orf(dna):
    """Searches all 6 frames for the longest sequence without a stop codon."""
    dna = dna.upper()
    strands = [dna, reverse_complement(dna)]
    all_orfs = []
    
    for strand in strands:
        for frame in range(3):
            translated = translate_dna(strand[frame:])
            orfs = translated.split('_')
            for o in orfs:
                if len(o) > 0:
                    all_orfs.append(o)
    
    if not all_orfs:
        return ""
        
    return max(all_orfs, key=len)

def prepare_fasta(input_csv, output_fasta, target_features=None):
    print("="*70)
    print("STEP 3: Structural Validation Prep (FASTA)")
    print(f"Input: {input_csv}")
    if target_features:
        print(f"Filter: Only processing Feature IDs: {target_features}")
    print("="*70)
    
    if not os.path.exists(input_csv):
        print(f"Error: {input_csv} not found.")
        return

    df = pd.read_csv(input_csv, encoding='utf-8-sig')
    
    # 1. Apply Target Feature Filter (Discovery Selection)
    if target_features:
        df = df[df['feature_id'].isin(target_features)].copy()
    
    # 2. Filter for high-confidence assembly
    if 'method' in df.columns:
        df = df[df['method'] == 'Consensus Assembly'].copy()

    if df.empty:
        print("No matching sequences found for these parameters.")
        return

    print(f"Processing {len(df)} selected features...")
    
    # Create FASTA package
    with open(output_fasta, "w") as f:
        count = 0
        for i, row in df.iterrows():
            dna = row['sequence']
            fid = row['feature_id']
            protein = find_longest_orf(dna)
            
            if len(protein) < 20: 
                print(f"   [SKIP] Feature {fid}: Protein too short ({len(protein)} aa)")
                continue
                
            f.write(f">feature_{fid}\n{protein}\n")
            print(f"   [READY] Feature {fid}: {len(protein)} amino acids")
            count += 1
    
    if count > 0:
        print(f"\nFASTA Created: {output_fasta}")
        print("NEXT STEP: Feed this file to AlphaFold 2 / 3 Colab or local cluster.")
    print("="*70)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PlatyGeno Step 3: 6-Frame FASTA Prep.")
    parser.add_argument("--input", type=str, default="PLG_Stage2_Validation.csv", help="Path to novelty CSV")
    parser.add_argument("--output", type=str, default="PLG_Stage3_Structural_Prep.fasta", help="Output FASTA path")
    parser.add_argument("--features", type=str, help="Optional: Comma-separated Feature IDs to extract (e.g. 7393,1254)")
    
    args = parser.parse_args()
    
    target_ids = [int(x.strip()) for x in args.features.split(",")] if args.features else None
    prepare_fasta(args.input, args.output, target_ids)
