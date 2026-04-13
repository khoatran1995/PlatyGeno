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
            # Split by STOP codons and find the longest segment
            orfs = translated.split('_')
            for o in orfs:
                if len(o) > 0:
                    all_orfs.append(o)
    
    if not all_orfs:
        return ""
        
    return max(all_orfs, key=len)

def prepare_fasta(input_csv, output_fasta):
    print("="*70)
    print("STEP 3: Turbo-ORF FASTA Preparation (6-Frame Detection)")
    print(f"Input: {input_csv}")
    print("="*70)
    
    if not os.path.exists(input_csv):
        print(f"Error: {input_csv} not found.")
        return

    df = pd.read_csv(input_csv, encoding='utf-8-sig')
    
    # Filter for high-confidence assembly if multiple exist
    if 'method' in df.columns:
        # We prefer Consensus Assembly for better AlphaFold results
        df = df[df['method'] == 'Consensus Assembly'].copy()

    print(f"Loaded {len(df)} candidate features.")

    # Create FASTA package
    with open(output_fasta, "w") as f:
        count = 0
        for i, row in df.iterrows():
            dna = row['sequence']
            fid = row['feature_id']
            
            # Find the best protein encoding
            protein = find_longest_orf(dna)
            
            # AlphaFold needs at least 20-30 amino acids for a meaningful fold
            if len(protein) < 20: 
                print(f"   [SKIP] Feature {fid}: Protein too short ({len(protein)} aa)")
                continue
                
            f.write(f">feature_{fid}\n{protein}\n")
            print(f"   [READY] Feature {fid}: {len(protein)} amino acids")
            count += 1
    
    if count > 0:
        print(f"\nFASTA package created: {output_fasta}")
        print("Ready for structural folding.")
    else:
        print("\nNo sequences were long enough for folding.")
        
    print("="*70)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PlatyGeno Step 3: 6-Frame FASTA Prep.")
    parser.add_argument("--input", type=str, required=True, help="Path to novelty CSV")
    parser.add_argument("--output", type=str, default="validation/potential_novel_sequences.fasta", help="Output FASTA path")
    
    args = parser.parse_args()
    prepare_fasta(args.input, args.output)
