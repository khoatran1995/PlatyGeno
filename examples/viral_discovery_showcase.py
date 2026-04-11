"""
PlatyGeno: Viral Discovery Showcase 🧬
--------------------------------------
This script demonstrates the end-to-end gene discovery pipeline:
1. Data Acquisition: Automatic download of a clinical virome sample.
2. Discovery: Unsupervised feature scanning via Evo 2 and SAEs.
3. Syntactic Validation: Automated NCBI BLAST search for novelty detection.
4. Structural Validation: 3D Protein Folding via Meta ESMFold API.

Developed by Khoa Tu Tran for PhD Technical Benchmarking.
"""

import os
import gzip
import shutil
import pandas as pd
import requests
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import platygeno

# --- 1. BIOLOGICAL UTILITIES ---

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
    """Submits DNA to NCBI BLAST. Returns (top_hit_title, identity_perc, is_novel)."""
    print(f"   🔎 Querying NCBI BLAST...")
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", dna_sequence, format_type="XML")
        blast_record = NCBIXML.read(result_handle)
        if len(blast_record.alignments) > 0:
            top_hit = blast_record.alignments[0]
            top_hsp = top_hit.hsps[0]
            identity = (top_hsp.identities / top_hsp.align_length) * 100
            is_novel = identity < 70
            title = top_hit.title[:50] + "..." if len(top_hit.title) > 50 else top_hit.title
            return title, identity, is_novel
        return "No Hits Found", 0, True
    except Exception as e:
        return f"BLAST Error: {str(e)[:20]}", 0, False

def fold_protein(sequence):
    """Calls Meta ESMFold API. Returns (pdb_string, mean_plddt)."""
    url = "https://api.esmatlas.com/fold/v1/pdb/"
    try:
        response = requests.post(url, data=sequence, timeout=90)
        pdb_string = response.text
        plddts = [float(line[60:66].strip()) for line in pdb_string.splitlines() if line.startswith("ATOM") and " CA " in line]
        mean_plddt = sum(plddts) / len(plddts) if plddts else 0
        return pdb_string, mean_plddt
    except Exception as e:
        print(f"      ⚠️ Folding API Error: {e}")
        return None, 0

# --- 2. DATA HANDLING ---

def download_and_extract():
    """Handles automatic data acquisition and validation."""
    raw_path = "data/raw/benchmark_virome.fastq.gz"
    fastq_path = "data/raw/benchmark_virome.fastq"
    url = "https://github.com/beard-group/virome-benchmarking/raw/master/data/small_virome.fastq.gz"
    
    if not os.path.exists(raw_path) or os.path.getsize(raw_path) < 100:
        print(f"📥 Downloading Clinical Virome Sample...")
        os.makedirs("data/raw", exist_ok=True)
        resp = requests.get(url, stream=True)
        with open(raw_path, 'wb') as f:
            shutil.copyfileobj(resp.raw, f)
            
    if not os.path.exists(fastq_path) or os.path.getsize(fastq_path) == 0:
        print(f"📦 Extracting Genomic Data...")
        with gzip.open(raw_path, 'rb') as f_in, open(fastq_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    return fastq_path

# --- 3. MASTER PIPELINE ---

def run_showcase():
    print("="*60)
    print("🧬 PLATYGENO VIRAL DISCOVERY SHOWCASE")
    print("="*60)
    
    fastq_path = download_and_extract()
    print(f"📊 Dataset Ready: {os.path.getsize(fastq_path)/1024:.1f} KB")

    # Start Discovery
    print("\n🔍 Phase 1: Unsupervised Feature Scanning...")
    results = platygeno.discover_genes(
        input_path=fastq_path,
        scan_end=None, # Full Scan
        min_activation=5.0,
        top_n=10
    )
    
    if results.empty:
        print("❌ No significant features discovered in this sample.")
        return

    print(f"\n✅ Phase 2: {len(results)} Features Isolated. Starting Validation...")
    
    final_report = []
    os.makedirs("data/folds", exist_ok=True)
    
    # Process the top 5 most unique discoveries
    for i, row in results.head(5).iterrows():
        dna_seq = row['sequence']
        feature_id = row['feature_id']
        
        print(f"\n➡️ Discovery #{i+1} [Feature ID: {feature_id}]")
        
        # Validation A: Syntactic (BLAST)
        hit, identity, is_novel = blast_dna(dna_seq)
        status = "🌟 NOVEL" if is_novel else "🧬 REDISCOVERY"
        print(f"   [{status}] BLAST: {hit} ({identity:.1f}%)")
        
        # Validation B: Structural (Folding)
        protein_seq = translate_dna(dna_seq)
        print(f"   🧬 Folding Protein ({len(protein_seq)}aa)...")
        pdb_data, plddt = fold_protein(protein_seq)
        
        if pdb_data:
            with open(f"data/folds/feature_{feature_id}.pdb", "w") as f:
                f.write(pdb_data)
            print(f"   ✨ Structural Confidence: {plddt:.2f}")
            
            final_report.append({
                'Rank': i+1,
                'FeatureID': feature_id,
                'Novelty': status,
                'BLAST_Hit': hit,
                'Identity': f"{identity:.1f}%",
                'pLDDT': f"{plddt:.2f}",
                'ProteinLen': len(protein_seq)
            })

    # --- 4. TECHNICAL REPORT GENERATION ---
    if final_report:
        report_df = pd.DataFrame(final_report)
        report_df.to_csv("data/benchmark_results_summary.csv", index=False)
        
        print("\n\n" + "="*60)
        print("🏆 FINAL TECHNICAL DISCOVERY REPORT")
        print("="*60)
        print(report_df.to_markdown(index=False))
        print("="*60)
        print("\nSummary: All PDB models are saved in data/folds/.")
        print("Copy the table above for your PhD Technical Report.")

if __name__ == "__main__":
    run_showcase()
