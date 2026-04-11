import os
import gzip
import shutil
import pandas as pd
import requests
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import platygeno

# --- Scientific Utilities ---

def translate_dna(dna):
    """Translates DNA to protein across all 6 frames and returns the longest ORF."""
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
    """Submits DNA to NCBI BLAST and returns novelty metrics."""
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
    """Folds protein using Meta ESMFold API and extracts pLDDT accuracy."""
    url = "https://api.esmatlas.com/fold/v1/pdb/"
    try:
        response = requests.post(url, data=sequence, timeout=60)
        pdb_string = response.text
        plddts = [float(line[60:66].strip()) for line in pdb_string.splitlines() if line.startswith("ATOM") and " CA " in line]
        return pdb_string, (sum(plddts) / len(plddts) if plddts else 0)
    except Exception:
        return None, 0

# --- Data Management ---

def get_discovery_dataset(dest_path):
    """Ensures a valid dataset exists at the destination."""
    url = "https://zenodo.org/record/582600/files/mutant_R1.fastq"
    if not os.path.exists(dest_path) or os.path.getsize(dest_path) == 0:
        print(f"📥 Loading discovery dataset from Zenodo ({os.path.basename(dest_path)})...")
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        try:
            response = requests.get(url, stream=True)
            with open(dest_path, 'wb') as f:
                shutil.copyfileobj(response.raw, f)
            print("✅ Dataset acquisition successful.")
        except Exception as e:
            print(f"❌ Failed to download dataset: {e}")
            return False
    return True

# --- Main Showcase ---

def run_discovery_showcase():
    print("🧬 Starting PlatyGeno PhD Discovery Benchmark...")
    input_path = "data/raw/benchmark_sample.fastq"
    
    if not get_discovery_dataset(input_path):
        return

    # 1. Scanning
    print(f"🔍 Analyzing dataset ({os.path.getsize(input_path)/1024:.1f} KB)...")
    results = platygeno.discover_genes(
        input_path=input_path,
        scan_end=None, # Full file scan
        min_activation=5.0,
        top_n=10
    )
    
    if results.empty:
        print("⚠️ No unique genomic features detected. Scan ended.")
        return

    # 2. Validation Loop
    print(f"✅ Isolated {len(results)} high-confidence features. Validating...")
    final_report = []
    os.makedirs("data/folds", exist_ok=True)
    
    for i, row in results.iterrows():
        dna = row['sequence']
        
        # BLAST Check
        hit_title, identity, is_novel = blast_validation(dna)
        status = "🌟 NOVEL" if is_novel else "🧬 KNOWN"
        
        # ESMFold Check
        protein = translate_dna(dna)
        pdb_data, plddt = esm_fold_protein(protein)
        
        print(f"   [{status}] Feature {row['feature_id']} | pLDDT: {plddt:.1f} | Top Hit: {hit_title}")

        if pdb_data:
            with open(f"data/folds/feature_{row['feature_id']}.pdb", "w") as f:
                f.write(pdb_data)
        
        final_report.append({
            'feature_id': row['feature_id'],
            'novelty': status,
            'plddt': plddt,
            'blast_identity': identity,
            'top_hit': hit_title,
            'sequence': dna
        })

    # 3. Final Technical Report
    report_df = pd.DataFrame(final_report)
    report_df.to_csv("data/benchmarking_full_report.csv", index=False)
    print("\n🏆 BENCHMARK COMPLETE")
    print(report_df[['feature_id', 'novelty', 'plddt', 'blast_identity']])
    print("\nFull data saved to data/benchmarking_full_report.csv and data/folds/.")

if __name__ == "__main__":
    run_discovery_showcase()
