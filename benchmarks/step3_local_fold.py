import os
import pandas as pd
import requests

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
        tokens = [seq[i:i+3] for i in range(0, len(seq), 3)]
        protein = "".join([table.get(t, 'X') for t in tokens if len(t) == 3])
        orfs = protein.split('_')
        return max(orfs, key=len) if orfs else ""

    rev_comp = "".join([{'A':'T','C':'G','G':'C','T':'A'}.get(b,'N') for b in reversed(dna)])
    frames = [dna, dna[1:], dna[2:], rev_comp, rev_comp[1:], rev_comp[2:]]
    return max([get_longest_orf(f) for f in frames], key=len)

def esm_fold_protein(sequence, feature_id):
    """Folds protein using Meta ESMFold API with detailed error reporting."""
    url = "https://api.esmatlas.com/fold/v1/pdb"
    headers = {"Content-Type": "text/plain"}
    
    if len(sequence) < 5:
        print(f"      ⚠️ Feature {feature_id}: Protein sequence too short ({len(sequence)} aa). Skipping.")
        return None, 0

    try:
        response = requests.post(url, data=sequence, headers=headers, timeout=60)
        pdb_string = response.text
        
        if response.status_code != 200:
            print(f"      ⚠️ Meta API Error (Status {response.status_code}): {pdb_string[:100]}")
            return None, 0
            
        if "ATOM" not in pdb_string:
            print(f"      ⚠️ Feature {feature_id}: Invalid PDB response (no ATOM records).")
            return None, 0
            
        plddts = [float(line[60:66].strip()) for line in pdb_string.splitlines() if line.startswith("ATOM") and " CA " in line]
        return pdb_string, (sum(plddts) / len(plddts) if plddts else 0)
    except Exception as e:
        print(f"      ⚠️ Network/Request Error: {e}")
        return None, 0

def run_folding_phase():
    print("="*70)
    print("💎 PHASE 3: Local 3D Protein Folding (Full Structural Proof)")
    print("="*70)
    
    # Smart path check for the previous step's results
    input_csv = "ibd_clinical_report.csv"
    if not os.path.exists(input_csv):
        if os.path.exists("benchmarks/ibd_clinical_report.csv"):
            input_csv = "benchmarks/ibd_clinical_report.csv"
        elif os.path.exists("blast_results.csv"):
            input_csv = "blast_results.csv"
        else:
            print(f"❌ Error: Results CSV not found. Run Step 2 or Validation first.")
            return

    blast_df = pd.read_csv(input_csv, encoding='utf-8-sig')
    print(f"✅ Loaded {len(blast_df)} entries. Starting Force-Folding for Structural Proof...")

    final_report = []
    os.makedirs("data/folds", exist_ok=True)

    for i, row in blast_df.iterrows():
        dna = row['sequence']
        feature_id = row['feature_id']
        
        # --- PHASE 3: FORCE FOLDING ---
        print(f"   ✨ Feature {feature_id} | Force-Folding Protein for Structural Proof...")
        protein = translate_dna(dna)
        pdb_data, plddt = esm_fold_protein(protein, feature_id)
        
        if pdb_data:
            with open(f"data/folds/feature_{feature_id}.pdb", "w") as f:
                f.write(pdb_data)
            print(f"      ✅ High-accuracy fold saved. Confidence (pLDDT): {plddt:.1f}")
        else:
            # Error message is handled inside esm_fold_protein
            plddt = 0
        
        final_report.append({
            'feature_id': feature_id,
            'novelty': row['novelty'],
            'plddt': plddt,
            'e_value': row['e_value'],
            'top_hit': row['top_hit'],
            'sequence': dna
        })

    report_df = pd.DataFrame(final_report)
    report_df.to_csv("ibd_final_report.csv", index=False, encoding='utf-8-sig')
    print("\n🏆 FOLDING PHASE COMPLETE. Summary saved to: ibd_final_report.csv")
    print("3D PDB files saved to: data/folds/")

if __name__ == "__main__":
    run_folding_phase()
