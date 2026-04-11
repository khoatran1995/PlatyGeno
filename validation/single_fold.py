import sys
import os
import torch
import esm

def translate_dna(dna):
    """Simplified DNA translation for single-sequence folding."""
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

def main():
    if len(sys.argv) < 2:
        print("Usage: python single_fold.py <DNA_SEQUENCE>")
        return

    dna_seq = sys.argv[1].upper()
    print(f"🧬 DNA Input: {dna_seq}")
    
    protein_seq = translate_dna(dna_seq)
    print(f"💎 Translated Protein: {protein_seq} ({len(protein_seq)} aa)")

    print("⏳ Loading ESMFold model to GPU...")
    model = esm.pretrained.esmfold_v1().eval().cuda()

    print("🚀 Folding...")
    with torch.no_grad():
        pdb_content = model.infer_pdb(protein_seq)
    
    output_name = "single_discovery.pdb"
    with open(output_name, "w") as f:
        f.write(pdb_content)

    print(f"✅ Structural discovery saved to: {output_name}")

if __name__ == "__main__":
    main()
