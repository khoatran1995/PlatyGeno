import os
import argparse
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

def run_benchmarking_validation(input_csv="discovery_hits.csv", output_prefix=None):
    """Phase 2: Standard Novelty Validation (Supports Chunked Inputs)"""
    
    # Generate output names based on input name if prefix isn't provided
    if output_prefix is None:
        base_name = os.path.splitext(input_csv)[0]
        output_csv = f"blast_results_{base_name.replace('discovery_hits_', '')}.csv" if "discovery_hits_" in base_name else "blast_results.csv"
        novel_csv = f"novel_sequences_{base_name.replace('discovery_hits_', '')}.csv" if "discovery_hits_" in base_name else "potential_novel_sequences.csv"
    else:
        output_csv = f"blast_results_{output_prefix}.csv"
        novel_csv = f"novel_sequences_{output_prefix}.csv"

    print("="*70)
    print(f"PHASE 2: Remote NCBI BLAST Validation")
    print(f"Input:  {input_csv}")
    print(f"Audit:  {output_csv}")
    print("="*70)

    if not os.path.exists(input_csv):
        print(f"❌ Error: {input_csv} not found.")
        print(f"Make sure you specified the correct chunk file (e.g., discovery_hits_5001_10000.csv).")
        return

    # 1. Load data
    df = pd.read_csv(input_csv)
    # Focusing on consensus sequences (highest quality)
    if 'method' in df.columns:
        consensus_df = df[df['method'] == 'Consensus Assembly'].copy()
    else:
        consensus_df = df.copy()
        
    consensus_df = consensus_df.drop_duplicates(subset=['feature_id'])

    if consensus_df.empty:
        print("No valid sequences found for validation.")
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
    parser = argparse.ArgumentParser(description="PlatyGeno Step 2: BLAST Validation.")
    parser.add_argument("--input", type=str, default="discovery_hits.csv", help="Discovery CSV to validate")
    parser.add_argument("--output-prefix", type=str, help="Prefix for result files")
    
    args = parser.parse_args()
    
    run_benchmarking_validation(input_csv=args.input, output_prefix=args.output_prefix)
