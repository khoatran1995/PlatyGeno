import os
import gzip
import shutil
import pandas as pd
import requests
import platygeno

def translate_dna(dna):
    """
    Translates DNA to protein across all 6 frames and returns the longest ORF.
    Simple but effective for the showcase.
    """
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
    
    def get_longest_orf_in_frame(seq):
        tokens = [seq[i:i+3] for i in range(0, len(seq), 3)]
        protein = "".join([table.get(t, 'X') for t in tokens if len(t) == 3])
        # Find longest orf between '_' or just the longest string
        orfs = protein.split('_')
        return max(orfs, key=len) if orfs else ""

    rev_comp = "".join([{'A':'T','C':'G','G':'C','T':'A'}.get(b,'N') for b in reversed(dna)])
    
    frames = [dna, dna[1:], dna[2:], rev_comp, rev_comp[1:], rev_comp[2:]]
    orfs = [get_longest_orf_in_frame(f) for f in frames]
    return max(orfs, key=len)

def fold_protein(sequence):
    """
    Calls Meta ESMFold API. Returns (pdb_string, mean_plddt).
    """
    url = "https://api.esmatlas.com/fold/v1/pdb/"
    try:
        response = requests.post(url, data=sequence, timeout=60)
        pdb_string = response.text
        
        # Extract pLDDT from B-factor column
        plddts = []
        for line in pdb_string.splitlines():
            if line.startswith("ATOM") and " CA " in line:
                plddt = float(line[60:66].strip())
                plddts.append(plddt)
        
        mean_plddt = sum(plddts) / len(plddts) if plddts else 0
        return pdb_string, mean_plddt
    except Exception as e:
        print(f"Error folding: {e}")
        return None, 0

def run_showcase():
    print("🧬 Starting PlatyGeno Viral Discovery Showcase...")
    
    raw_path = "data/raw/benchmark_virome.fastq.gz"
    fastq_path = "data/raw/benchmark_virome.fastq"
    
    # 1. Extraction
    if not os.path.exists(fastq_path):
        print("📦 Extracting sample dataset...")
        with gzip.open(raw_path, 'rb') as f_in:
            with open(fastq_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    # 2. Discovery
    print("🔍 Scanning for novel viral features...")
    # NOTE: In a real run, threshold 8.0-10.0 is best for "rare" signals
    results = platygeno.discover_genes(
        input_path=fastq_path,
        scan_end=2000,
        min_activation=8.0,
        top_n=5
    )
    
    if results.empty:
        print("⚠️ No strong signals found in this subsample. Try lowering threshold or increasing scan_end.")
        return

    # 3. Validation & Folding
    print(f"✅ Found {len(results)} potential features. Folding top 3...")
    
    final_report = []
    os.makedirs("data/folds", exist_ok=True)
    
    for i, row in results.head(3).iterrows():
        dna_seq = row['sequence']
        protein_seq = translate_dna(dna_seq)
        
        print(f"➡️ Folding Feature {row['feature_id']} (Protein Len: {len(protein_seq)}aa)...")
        pdb_data, plddt = fold_protein(protein_seq)
        
        if pdb_data:
            pdb_filename = f"data/folds/feature_{row['feature_id']}.pdb"
            with open(pdb_filename, "w") as f:
                f.write(pdb_data)
            
            final_report.append({
                'feature_id': row['feature_id'],
                'dna_len': len(dna_seq),
                'p_len': len(protein_seq),
                'plddt': plddt,
                'pdb_file': pdb_filename
            })
            print(f"   🌟 Confidence (pLDDT): {plddt:.2f}")

    # 4. Save Final Report
    report_df = pd.DataFrame(final_report)
    report_df.to_csv("data/benchmark_results_summary.csv", index=False)
    print("\n🏆 SHOWCASE COMPLETE")
    print("="*40)
    print(report_df[['feature_id', 'plddt', 'p_len']])
    print("="*40)
    print("All PDB structures and the summary report are in data/folds and data/benchmark_results_summary.csv.")

if __name__ == "__main__":
    run_showcase()
