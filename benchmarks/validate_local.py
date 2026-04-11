import os
import pandas as pd
import requests
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# --- Scientific Utilities (Local Version) ---

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

def blast_validation(dna_sequence):
    """Submits DNA to NCBI BLAST."""
    print(f"   🔎 Querying NCBI BLAST...")
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", dna_sequence, format_type="XML")
        blast_record = NCBIXML.read(result_handle)
        if len(blast_record.alignments) > 0:
            top_hit = blast_record.alignments[0]
            identity = (top_hit.hsps[0].identities / top_hit.hsps[0].align_length) * 100
            return top_hit.title[:50], identity, (identity < 70)
        return "No Hits Found", 0, True
    except Exception as e:
        return f"BLAST Error ({e})", 0, False

def esm_fold_protein(sequence):
    """Folds protein using Meta ESMFold API."""
    url = "https://api.esmatlas.com/fold/v1/pdb"
    headers = {"Content-Type": "text/plain"}
    try:
        response = requests.post(url, data=sequence, headers=headers, timeout=60)
        pdb_string = response.text
        if "ATOM" not in pdb_string:
            return None, 0
        plddts = [float(line[60:66].strip()) for line in pdb_string.splitlines() if line.startswith("ATOM") and " CA " in line]
        return pdb_string, (sum(plddts) / len(plddts) if plddts else 0)
    except Exception:
        return None, 0

def run_local_validation():
    print("="*70)
    print("🔬 Starting PlatyGeno Local Validation Phase (Free CPU Logic)")
    print("="*70)
    
    # Smart path check: looks in root, then in the benchmarks/ folder
    input_csv = "discovery_hits.csv"
    if not os.path.exists(input_csv):
        if os.path.exists("benchmarks/discovery_hits.csv"):
            input_csv = "benchmarks/discovery_hits.csv"
        else:
            print(f"❌ Error: {input_csv} not found in root or benchmarks/ folder.")
            return

    hits_df = pd.read_csv(input_csv)
    print(f"✅ Loaded {len(hits_df)} unique motifs from RunPod. Starting Scientific Validation...")

    final_report = []
    os.makedirs("data/folds", exist_ok=True)

    for i, row in hits_df.iterrows():
        dna = row['sequence']
        feature_id = row['feature_id']
        
        # 1. BLAST (Run for all 100)
        hit_title, identity, is_novel = blast_validation(dna)
        status = "🌟 NOVEL" if is_novel else "🧬 KNOWN"
        
        # 2. Conditional ESMFold (Only for NEW Genes)
        pdb_data, plddt = None, 0
        if is_novel:
            print(f"   [{status}] Feature {feature_id} | Folding Protein (Discovery Mode)...")
            protein = translate_dna(dna)
            pdb_data, plddt = esm_fold_protein(protein)
            
            if pdb_data:
                with open(f"data/folds/feature_{feature_id}.pdb", "w") as f:
                    f.write(pdb_data)
                print(f"   ✨ Confidence Score (pLDDT): {plddt:.2f}")
        else:
            print(f"   [{status}] Feature {feature_id} | Skipping ESMFold (Rediscovery). Hit: {hit_title}")
        
        final_report.append({
            'feature_id': feature_id,
            'discovery_type': row.get('discovery_type', 'N/A'),
            'method': row.get('method', 'N/A'),
            'novelty': status,
            'plddt': plddt,
            'blast_identity': identity,
            'top_hit': hit_title,
            'sequence': dna
        })

    # 3. Final Report
    report_df = pd.DataFrame(final_report)
    report_df.to_csv("ibd_clinical_report.csv", index=False)
    print("\n🏆 VALIDATION COMPLETE")
    print("="*70)
    print(report_df[['feature_id', 'novelty', 'plddt', 'blast_identity']])
    print("="*70)
    print("\nFull PhD Technical Report saved as: ibd_clinical_report.csv")

if __name__ == "__main__":
    run_local_validation()
