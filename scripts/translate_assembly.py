import sys

def translate(dna):
    dna = dna.upper().replace('U', 'T')
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
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        protein += table.get(codon, 'X')
    return protein

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python translate_assembly.py <fasta_file>")
        sys.exit(1)
    
    with open(sys.argv[1], 'r') as f:
        lines = f.readlines()
        dna = "".join([l.strip() for l in lines if not l.startswith('>')])
    
    # Try all 3 forward frames
    print(f"Translating {sys.argv[1]}...")
    for frame in range(3):
        prot = translate(dna[frame:])
        print(f"\nFrame {frame+1}:")
        print(prot)
