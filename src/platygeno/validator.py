import os
import pandas as pd
import time
import socket
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# Set a global timeout for socket operations to prevent indefinite hanging
socket.setdefaulttimeout(300) # 5 minutes

def blast_query(dna_sequence, max_retries=3):
    """
    Submits DNA to NCBI BLAST with retry logic.
    Low-level utility for the validator.
    """
    for attempt in range(max_retries):
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
            time.sleep(30)
    return "BLAST Failed", 0, False, 1.0

def validate_novelty(df, output_path=None, novel_path=None):
    """
    Library-native validator for PlatyGeno discovery results.
    Processes a Discovery DataFrame and identifies 'Novel' vs 'Known' sequences.
    
    Args:
        df (pd.DataFrame): The discovery hits from discover_genes().
        output_path (str, optional): Path to save the full audit trail.
        novel_path (str, optional): Path to save potential novel sequences.
        
    Returns:
        pd.DataFrame: A DataFrame enriched with BLAST validation results.
    """
    print("="*70)
    print("PlatyGeno Validator: Local NCBI BLAST Novelty Check")
    print("="*70)
    
    if df.empty:
        print("⚠️ Input DataFrame is empty. Nothing to validate.")
        return df

    # 1. Check for existing progress if output_path is provided (Checkpointing)
    processed_ids = set()
    if output_path and os.path.exists(output_path):
        try:
            existing_df = pd.read_csv(output_path, encoding='utf-8-sig')
            if 'feature_id' in existing_df.columns:
                processed_ids = set(existing_df['feature_id'].unique())
                print(f"   Found existing progress. Resuming from {len(processed_ids)} hits.")
        except Exception:
            pass

    # 2. Filtering Strategy: Focus on Consensus Assemblies
    if 'method' in df.columns:
        consensus_df = df[df['method'] == 'Consensus Assembly'].copy()
        if consensus_df.empty:
            consensus_df = df.copy()
    else:
        consensus_df = df.copy()
        
    # Deduplicate and Filter
    consensus_df = consensus_df.drop_duplicates(subset=['feature_id'])
    remaining_df = consensus_df[~consensus_df['feature_id'].isin(processed_ids)]
    
    if remaining_df.empty:
        print("   All features have already been validated.")
        return df

    print(f"   Starting BLAST validation for {len(remaining_df)} unique consensus sequences...")

    # 3. Validation Loop
    for i, row in remaining_df.iterrows():
        dna = row['sequence']
        feature_id = row['feature_id']
        
        hit_title, identity, is_novel, e_value = blast_query(dna)
        status = "NOVEL" if is_novel else "KNOWN"
        
        print(f"      [{status}] Feature {feature_id} | Identity: {identity:.1f}% | E-val: {e_value}")

        # Update the main results batch
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
        
        # Incremental Save (Checkpointing)
        if output_path:
            file_exists = os.path.exists(output_path)
            res_df.to_csv(output_path, mode='a', index=False, header=not file_exists, encoding='utf-8-sig')

    # 4. Final Processing & Summary
    if output_path:
        final_df = pd.read_csv(output_path, encoding='utf-8-sig')
        if novel_path:
            novel_df = final_df[final_df['novelty'] == "NOVEL"]
            novel_df.to_csv(novel_path, index=False, encoding='utf-8-sig')
            print(f"💾 Novel sequences saved to: {novel_path}")
        return final_df
    
    return pd.DataFrame() # If no path provided, just return empty (path is required for this version)
