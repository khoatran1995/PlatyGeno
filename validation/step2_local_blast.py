import os
import pandas as pd
import time
import socket
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# Set a global timeout for socket operations to prevent indefinite hanging
socket.setdefaulttimeout(300) # 5 minutes

def blast_validation(dna_sequence, max_retries=3):
    """Submits DNA to NCBI BLAST with retry logic."""
    for attempt in range(max_retries):
        print(f"   Querying NCBI BLAST (Attempt {attempt + 1}/{max_retries})...")
        try:
            if attempt > 0:
                time.sleep(30) # Wait longer between retries
                
            result_handle = NCBIWWW.qblast("blastn", "nt", dna_sequence, format_type="XML")
            blast_record = NCBIXML.read(result_handle)
            
            if len(blast_record.alignments) > 0:
                top_hit = blast_record.alignments[0]
                hsp = top_hit.hsps[0]
                identity = (hsp.identities / hsp.align_length) * 100
                e_value = hsp.expect
                time.sleep(15) # Polite delay for NCBI servers
                return top_hit.title[:50], identity, (identity < 70), e_value
            return "No Hits Found", 0, True, 10.0
        except Exception as e:
            print(f"      BLAST error (attempt {attempt+1}): {e}")
            time.sleep(30)
    return "BLAST Failed", 0, False, 1.0

def regenerate_novel_list(full_path, novel_path):
    """Refreshes the potential_novel_sequences.csv from the main audit trail."""
    if os.path.exists(full_path):
        try:
            full_df = pd.read_csv(full_path, encoding='utf-8-sig')
            novel_df = full_df[full_df['novelty'] == "NOVEL"]
            novel_df.to_csv(novel_path, index=False, encoding='utf-8-sig')
        except Exception as e:
            print(f"Warning: Could not regenerate novel list: {e}")

def run_blast_phase():
    print("="*70)
    print("PHASE 2: Local NCBI BLAST Validation (Checkpoint-Enabled)")
    print("="*70)
    
    # Path configuration
    base_dir = os.path.dirname(__file__)
    input_csv = "discovery_hits.csv"
    if not os.path.exists(input_csv):
        if os.path.exists("validation/discovery_hits.csv"):
            input_csv = "validation/discovery_hits.csv"
        else:
            print(f"Error: {input_csv} not found.")
            return

    full_audit_path = os.path.join(base_dir, "blast_results.csv")
    novel_path = os.path.join(base_dir, "potential_novel_sequences.csv")

    hits_df = pd.read_csv(input_csv, encoding='utf-8-sig')
    
    # 1. Load existing progress to resume
    processed_ids = set()
    if os.path.exists(full_audit_path):
        try:
            existing_df = pd.read_csv(full_audit_path, encoding='utf-8-sig')
            if 'feature_id' in existing_df.columns:
                processed_ids = set(existing_df['feature_id'].unique())
                print(f"Found existing progress. Skipping {len(processed_ids)} already processed features.")
        except Exception:
            print("Could not read existing audit trail. Starting fresh.")

    # 2. Filter for Consensus Assembly
    if 'method' in hits_df.columns:
        consensus_df = hits_df[hits_df['method'] == 'Consensus Assembly'].copy()
        if consensus_df.empty:
            consensus_df = hits_df.copy()
    else:
        consensus_df = hits_df.copy()
        
    # Deduplicate and filter out already processed IDs
    consensus_df = consensus_df.drop_duplicates(subset=['feature_id'])
    remaining_df = consensus_df[~consensus_df['feature_id'].isin(processed_ids)]
    
    if remaining_df.empty:
        print("All unique clinical features have already been processed.")
        regenerate_novel_list(full_audit_path, novel_path)
        return

    print(f"Starting BLAST searches for {len(remaining_df)} remaining features...")

    # 3. Process incrementally
    for i, row in remaining_df.iterrows():
        dna = row['sequence']
        feature_id = row['feature_id']
        
        hit_title, identity, is_novel, e_value = blast_validation(dna)
        status = "NOVEL" if is_novel else "KNOWN"
        
        print(f"   [{status}] Feature {feature_id} | Identity: {identity:.1f}% | E-value: {e_value}")

        # Create a single-row DataFrame for incremental saving
        new_result = {
            'feature_id': feature_id,
            'discovery_type': row.get('discovery_type', 'N/A'),
            'novelty': status,
            'e_value': e_value,
            'blast_identity': identity,
            'top_hit': hit_title,
            'sequence': dna
        }
        res_df = pd.DataFrame([new_result])
        
        # Append to CSV immediately
        file_exists = os.path.exists(full_audit_path)
        res_df.to_csv(full_audit_path, mode='a', index=False, header=not file_exists, encoding='utf-8-sig')

    # 4. Final summary refresh
    regenerate_novel_list(full_audit_path, novel_path)
    
    print("\nSTEP 2 COMPLETE")
    print("="*70)
    print(f"Full Audit Trail: {full_audit_path}")
    print(f"Filtered Novel Hits: {novel_path}")
    print("NEXT STEP: Run 'python validation/step3_fasta_prep.py'")
    print("="*70)

if __name__ == "__main__":
    run_blast_phase()
