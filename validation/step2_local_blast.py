import os
import pandas as pd
import time
import socket
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# Set a global timeout for socket operations
socket.setdefaulttimeout(300) # 5 minutes

def blast_validation(dna_sequence, max_retries=3):
    """Submits DNA to NCBI BLAST with retry logic."""
    for attempt in range(max_retries):
        try:
            if attempt > 0:
                time.sleep(30)
                
            result_handle = NCBIWWW.qblast("blastn", "nt", dna_sequence, format_type="XML")
            blast_record = NCBIXML.read(result_handle)
            
            if len(blast_record.alignments) > 0:
                top_hit = blast_record.alignments[0]
                hsp = top_hit.hsps[0]
                identity = (hsp.identities / hsp.align_length) * 100
                e_value = hsp.expect
                time.sleep(15) 
                return top_hit.title[:50], identity, (identity < 70), e_value
            return "No Hits Found", 0, True, 10.0
        except Exception as e:
            time.sleep(30)
    return "BLAST Failed", 0, False, 1.0

def run_benchmarking_validation():
    """Phase 2: Standard Novelty Validation (Full Re-run Mode)"""
    input_csv = "discovery_hits.csv"
    output_csv = "blast_results.csv"
    novel_csv = "potential_novel_sequences.csv"

    print("="*70)
    print("PHASE 2: Remote NCBI BLAST Validation (Full Re-run Mode)")
    print("="*70)

    if not os.path.exists(input_csv):
        print(f"Error: {input_csv} not found. Run Step 1 first.")
        return

    # 1. Load data
    df = pd.read_csv(input_csv)
    # Focusing on consensus sequences (highest quality)
    consensus_df = df[df['method'] == 'Consensus Assembly'].copy()
    consensus_df = consensus_df.drop_duplicates(subset=['feature_id'])

    if consensus_df.empty:
        print("No consensus assemblies found in discovery_hits.csv.")
        return

    print(f"Starting fresh BLAST searches for {len(consensus_df)} features...")

    # 2. Validation Loop
    results_list = []
    for i, row in consensus_df.iterrows():
        dna = row['sequence']
        fid = row['feature_id']
        
        hit_title, identity, is_novel, e_value = blast_validation(dna)
        status = "NOVEL" if is_novel else "KNOWN"
        
        print(f"   [{status}] Feature {fid} | Identity: {identity:.1f}% | E-value: {e_value}")

        # Add to list
        results_list.append({
            'feature_id': fid,
            'novelty': status,
            'blast_identity': identity,
            'top_hit': hit_title,
            'sequence': dna,
            'e_value': e_value
        })

    # 3. Final reporting (Always Overwrites)
    final_df = pd.DataFrame(results_list)
    final_df.to_csv(output_csv, index=False, encoding='utf-8-sig')
    
    novel_df = final_df[final_df['novelty'] == "NOVEL"]
    novel_df.to_csv(novel_csv, index=False, encoding='utf-8-sig')

    print("\nSTEP 2 COMPLETE")
    print("="*70)
    print(f"Audit Trail: {output_csv}")
    print(f"Novel Hits:  {novel_csv}")
    print("NEXT STEP: Run 'python validation/step3_fasta_prep.py'")
    print("="*70)

if __name__ == "__main__":
    run_benchmarking_validation()
