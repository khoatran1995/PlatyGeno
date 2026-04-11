import os
import pandas as pd
import time
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def blast_validation(dna_sequence, max_retries=3):
    """Submits DNA to NCBI BLAST with retry logic."""
    for attempt in range(max_retries):
        print(f"   🔎 Querying NCBI BLAST (Attempt {attempt + 1}/{max_retries})...")
        try:
            if attempt > 0:
                time.sleep(20)
                
            result_handle = NCBIWWW.qblast("blastn", "nt", dna_sequence, format_type="XML")
            blast_record = NCBIXML.read(result_handle)
            
            if len(blast_record.alignments) > 0:
                top_hit = blast_record.alignments[0]
                hsp = top_hit.hsps[0]
                identity = (hsp.identities / hsp.align_length) * 100
                e_value = hsp.expect
                time.sleep(10) # Polite delay
                return top_hit.title[:50], identity, (identity < 70), e_value
            return "No Hits Found", 0, True, 10.0
        except Exception as e:
            print(f"      ⚠️ BLAST error: {e}")
            time.sleep(20)
    return "BLAST Failed", 0, False, 1.0

def run_blast_phase():
    print("="*70)
    print("🔍 PHASE 2: Local NCBI BLAST Validation")
    print("="*70)
    
    input_csv = "discovery_hits.csv"
    if not os.path.exists(input_csv):
        if os.path.exists("benchmarks/discovery_hits.csv"):
            input_csv = "benchmarks/discovery_hits.csv"
        else:
            print(f"❌ Error: {input_csv} not found.")
            return

    hits_df = pd.read_csv(input_csv)
    print(f"✅ Loaded {len(hits_df)} hits. Starting BLAST searches...")

    blast_results = []
    for i, row in hits_df.iterrows():
        dna = row['sequence']
        feature_id = row['feature_id']
        
        hit_title, identity, is_novel, e_value = blast_validation(dna)
        status = "🌟 NOVEL" if is_novel else "🧬 KNOWN"
        
        print(f"   [{status}] Feature {feature_id} | Identity: {identity:.1f}% | E-value: {e_value}")

        blast_results.append({
            'feature_id': feature_id,
            'discovery_type': row.get('discovery_type', 'N/A'),
            'novelty': status,
            'e_value': e_value,
            'blast_identity': identity,
            'top_hit': hit_title,
            'sequence': dna
        })

    out_df = pd.DataFrame(blast_results)
    out_df.to_csv("blast_results.csv", index=False)
    print("\n🏆 BLAST PHASE COMPLETE. Results saved to: blast_results.csv")

if __name__ == "__main__":
    run_blast_phase()
