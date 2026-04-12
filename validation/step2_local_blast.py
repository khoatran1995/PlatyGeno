import os
import pandas as pd
import time
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def blast_validation(dna_sequence, max_retries=3):
    """Submits DNA to NCBI BLAST with retry logic."""
    for attempt in range(max_retries):
        print(f"   Querying NCBI BLAST (Attempt {attempt + 1}/{max_retries})...")
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
            print(f"      BLAST error: {e}")
            time.sleep(20)
    return "BLAST Failed", 0, False, 1.0

def run_blast_phase():
    print("="*70)
    print("PHASE 2: Local NCBI BLAST Validation")
    print("="*70)
    
    input_csv = "discovery_hits.csv"
    if not os.path.exists(input_csv):
        if os.path.exists("validation/discovery_hits.csv"):
            input_csv = "validation/discovery_hits.csv"
        else:
            print(f"Error: {input_csv} not found.")
            return

    hits_df = pd.read_csv(input_csv, encoding='utf-8-sig')
    
    # OPTIMIZATION: Only BLAST 'Consensus Assembly' sequences.
    # These have higher info density and better BLAST specificity than 30bp snippets.
    if 'method' in hits_df.columns:
        consensus_df = hits_df[hits_df['method'] == 'Consensus Assembly'].copy()
        # Fallback if no consensus found (should not happen in this pipeline)
        if consensus_df.empty:
            consensus_df = hits_df.copy()
    else:
        consensus_df = hits_df.copy()
        
    # Deduplicate by feature_id to avoid redundant queries
    consensus_df = consensus_df.drop_duplicates(subset=['feature_id'])
    
    print(f"Filtered to {len(consensus_df)} unique Consensus sequences. Starting BLAST searches...")

    blast_results = []
    for i, row in consensus_df.iterrows():
        dna = row['sequence']
        feature_id = row['feature_id']
        
        hit_title, identity, is_novel, e_value = blast_validation(dna)
        status = "NOVEL" if is_novel else "KNOWN"
        
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
    
    # Path standardization: Always save to the validation folder
    base_dir = os.path.dirname(__file__)
    full_path = os.path.join(base_dir, "blast_results.csv")
    novel_path = os.path.join(base_dir, "potential_novel_sequences.csv")
    
    # 1. Save Full Audit Trail
    out_df.to_csv(full_path, index=False, encoding='utf-8-sig')
    
    # 2. Filter and Save Potential Novel Sequences (Identity < 70% or no hits)
    novel_df = out_df[out_df['novelty'] == "NOVEL"]
    novel_df.to_csv(novel_path, index=False, encoding='utf-8-sig')
    
    print("\nSTEP 2 COMPLETE")
    print("="*70)
    print(f"Full Results: {full_path}")
    print(f"Novel Hits: {novel_path}")
    print(f"NEXT STEP: Run 'python validation/step3_fasta_prep.py'")
    print("="*70)

if __name__ == "__main__":
    run_blast_phase()
