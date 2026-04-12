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
    """Phase 2: Checkpoint-Enabled Novelty Validation"""
    input_csv = "discovery_hits.csv"
    output_csv = "blast_results.csv"
    novel_csv = "potential_novel_sequences.csv"

    print("="*70)
    print("PHASE 2: Remote NCBI BLAST Validation (Checkpoint-Enabled)")
    print("="*70)

    if not os.path.exists(input_csv):
        print(f"Error: {input_csv} not found. Run Step 1 first.")
        return

    # 1. Load data and identify pending work
    df = pd.read_csv(input_csv)
    # Focusing on consensus sequences (highest quality)
    consensus_df = df[df['method'] == 'Consensus Assembly'].copy()
    consensus_df = consensus_df.drop_duplicates(subset=['feature_id'])

    processed_ids = set()
    if os.path.exists(output_csv):
        try:
            existing_df = pd.read_csv(output_csv, encoding='utf-8-sig')
            processed_ids = set(existing_df['feature_id'].unique())
            print(f"Found existing progress. Resuming from {len(processed_ids)} features.")
        except Exception:
            pass

    remaining_df = consensus_df[~consensus_df['feature_id'].isin(processed_ids)]
    if remaining_df.empty:
        print("All features have already been validated.")
        return

    print(f"Starting BLAST searches for {len(remaining_df)} remaining features...")

    # 2. Validation Loop with Incremental Saving
    for i, row in remaining_df.iterrows():
        dna = row['sequence']
        fid = row['feature_id']
        
        hit_title, identity, is_novel, e_value = blast_validation(dna)
        status = "NOVEL" if is_novel else "KNOWN"
        
        print(f"   [{status}] Feature {fid} | Identity: {identity:.1f}% | E-value: {e_value}")

        # Save this result immediately
        result_row = {
            'feature_id': fid,
            'novelty': status,
            'blast_identity': identity,
            'top_hit': hit_title,
            'sequence': dna,
            'e_value': e_value
        }
        res_df = pd.DataFrame([result_row])
        file_exists = os.path.exists(output_csv)
        res_df.to_csv(output_csv, mode='a', index=False, header=not file_exists, encoding='utf-8-sig')

    # 3. Final reporting
    final_df = pd.read_csv(output_csv, encoding='utf-8-sig')
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
