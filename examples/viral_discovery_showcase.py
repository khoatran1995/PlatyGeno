import os
import gzip
import shutil
import pandas as pd
import requests
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import platygeno

def translate_dna(dna):
    """
    Translates DNA to protein across all 6 frames and returns the longest ORF.
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
        orfs = protein.split('_')
        return max(orfs, key=len) if orfs else ""

    rev_comp = "".join([{'A':'T','C':'G','G':'C','T':'A'}.get(b,'N') for b in reversed(dna)])
    frames = [dna, dna[1:], dna[2:], rev_comp, rev_comp[1:], rev_comp[2:]]
    orfs = [get_longest_orf_in_frame(f) for f in frames]
    return max(orfs, key=len)

def blast_dna(dna_sequence):
    """
    Submits DNA to NCBI BLAST. Returns (top_hit_title, identity_perc, is_novel).
    """
    print(f"   🔎 Querying NCBI BLAST (this may take 1-2 mins)...")
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", dna_sequence, format_type="XML")
        blast_record = NCBIXML.read(result_handle)
        if len(blast_record.alignments) > 0:
            top_hit = blast_record.alignments[0]
            top_hsp = top_hit.hsps[0]
            identity = (top_hsp.identities / top_hsp.align_length) * 100
            is_novel = identity < 70
            return top_hit.title[:50], identity, is_novel
        return "No Hits Found", 0, True
    except Exception as e:
        print(f"   ⚠️ BLAST Error: {e}")
        return "BLAST Timeout/Error", 0, False

def fold_protein(sequence):
    """
    Calls Meta ESMFold API. Returns (pdb_string, mean_plddt).
    """
    url = "https://api.esmatlas.com/fold/v1/pdb/"
    try:
        response = requests.post(url, data=sequence, timeout=60)
        pdb_string = response.text
        plddts = [float(line[60:66].strip()) for line in pdb_string.splitlines() if line.startswith("ATOM") and " CA " in line]
        mean_plddt = sum(plddts) / len(plddts) if plddts else 0
        return pdb_string, mean_plddt
    except Exception as e:
        print(f"Error folding: {e}")
        return None, 0

def is_valid_gzip(path):
    if not os.path.exists(path) or os.path.getsize(path) < 100:
        return False
    try:
        with gzip.open(path, 'rb') as f:
            f.read(10)
        return True
    except:
        return False

def download_dataset(dest_path):
    url = "https://github.com/beard-group/virome-benchmarking/raw/master/data/small_virome.fastq.gz"
    if not is_valid_gzip(dest_path):
        if os.path.exists(dest_path):
            print("⚠️ Data file corrupted. Re-downloading...")
            os.remove(dest_path)
        print(f"📥 Downloading benchmark dataset from GitHub...")
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        response = requests.get(url, stream=True)
        with open(dest_path, 'wb') as f:
            shutil.copyfileobj(response.raw, f)
        print("✅ Download complete.")

def run_showcase():
    print("🧬 Starting PlatyGeno Viral Discovery Showcase...")
    raw_path = "data/raw/benchmark_virome.fastq.gz"
    fastq_path = "data/raw/benchmark_virome.fastq"
    
    download_dataset(raw_path)
    if not os.path.exists(fastq_path):
        print("📦 Extracting sample dataset...")
        with gzip.open(raw_path, 'rb') as f_in, open(fastq_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    print("🔍 Scanning for novel viral features...")
    
    # Check file size for diagnostics
    file_size = os.path.getsize(fastq_path)
    print(f"📊 Dataset Size: {file_size/1024:.2f} KB")

    # Lower sensitivity to 5.0 and increase range to 4000 for better discovery
    results = platygeno.discover_genes(
        input_path=fastq_path, 
        scan_end=4000, 
        min_activation=5.0, 
        top_n=10
    )
    
    if results.empty:
        print("⚠️ No strong signals found. Try increasing scan_end.")
        return

    print(f"✅ Found {len(results)} potential features. Starting Parallel Validation...")
    final_report = []
    os.makedirs("data/folds", exist_ok=True)
    
    for i, row in results.iterrows():
        dna_seq = row['sequence']
        hit_title, identity, is_novel = blast_dna(dna_seq)
        protein_seq = translate_dna(dna_seq)
        print(f"   [{'🌟 NOVEL' if is_novel else '🧬 KNOWN'}] Feature {row['feature_id']} | plddt check...")
        pdb_data, plddt = fold_protein(protein_seq)
        
        if pdb_data:
            with open(f"data/folds/feature_{row['feature_id']}.pdb", "w") as f:
                f.write(pdb_data)
            final_report.append({'feature_id': row['feature_id'], 'dna_len': len(dna_seq), 'blast_hit': hit_title, 'identity': identity, 'is_novel': is_novel, 'plddt': plddt})
            print(f"   ✨ Confidence: {plddt:.2f}")

    report_df = pd.DataFrame(final_report)
    report_df.to_csv("data/benchmark_results_summary.csv", index=False)
    print("\n🏆 SHOWCASE COMPLETE")
    print(report_df[['feature_id', 'is_novel', 'identity', 'plddt']])

if __name__ == "__main__":
    run_showcase()
