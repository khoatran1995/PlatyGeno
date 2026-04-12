import os
import argparse
import pandas as pd
import time
import socket
import concurrent.futures
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# Set a global timeout for socket operations
socket.setdefaulttimeout(300) # 5 minutes

def blast_validation(dna_sequence, max_retries=3):
    """Submits DNA to NCBI BLAST with retry logic."""
    # We add a small staggered delay per thread to be polite to NCBI
    time.sleep(2) 
    
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
                return top_hit.title[:50], identity, (identity < 70), e_value
            return "No Hits Found", 0, True, 10.0
        except Exception as e:
            time.sleep(30)
    return "BLAST Failed", 0, False, 1.0

def process_feature_blast(row):
    """Worker function for parallel BLAST processing."""
    dna = row['sequence']
    fid = row['feature_id']
    
    hit_title, identity, is_novel, e_value = blast_validation(dna)
    status = "NOVEL" if is_novel else "KNOWN"
    
    print(f"   [{status}] Feature {fid} | Identity: {identity:.1f}% | E-value: {e_value}")
    
    return {
        'feature_id': fid,
        'novelty': status,
        'blast_identity': identity,
        'top_hit': hit_title,
        'sequence': dna,
        'e_value': e_value
    }

def run_benchmarking_validation(input_csv="discovery_hits.csv", output_prefix=None, max_workers=5, validate_all=False):
    """Phase 2: Turbo-BLAST Validation (Multi-threaded)"""
    
    # Generate output names
    if output_prefix is None:
        base_name = os.path.splitext(os.path.basename(input_csv))[0]
        suffix = base_name.replace('discovery_hits_', '')
        output_csv = f"blast_results_{suffix}.csv"
        novel_csv = f"novel_sequences_{suffix}.csv"
    else:
        output_csv = f"blast_results_{output_prefix}.csv"
        novel_csv = f"novel_sequences_{output_prefix}.csv"

    print("="*70)
    print(f"PHASE 2: Turbo-BLAST Validation (Threads: {max_workers})")
    print(f"Input:  {input_csv}")
    print(f"Mode:   {'Validate ALL (38+)' if validate_all else 'Consensus-Only (13)'}")
    print(f"Audit:  {output_csv}")
    print("="*70)

    if not os.path.exists(input_csv):
        print(f"❌ Error: {input_csv} not found.")
        return

    # 1. Load data
    df = pd.read_csv(input_csv)
    if not validate_all and 'method' in df.columns:
        valid_df = df[df['method'] == 'Consensus Assembly'].copy()
    else:
        valid_df = df.copy()
        
    if valid_df.empty:
        print("No valid sequences found for validation.")
        return

    print(f"Starting Turbo-BLAST for {len(valid_df)} features...")

    # 2. Parallel Validation Loop
    results_list = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Convert df to list of dicts for executor
        rows = valid_df.to_dict('records')
        future_to_info = {executor.submit(process_feature_blast, row): row for row in rows}
        
        for future in concurrent.futures.as_completed(future_to_info):
            row_orig = future_to_info[future]
            try:
                res = future.result()
                # Include the method metadata
                res['method'] = row_orig.get('method', 'Unknown')
                results_list.append(res)
            except Exception as exc:
                fid = row_orig.get('feature_id', '?')
                print(f"   [ERROR] Feature {fid} generated an exception: {exc}")

    # 3. Final reporting
    if not results_list:
        print("❌ All BLAST searches failed.")
        return

    final_df = pd.DataFrame(results_list)
    final_df.to_csv(output_csv, index=False, encoding='utf-8-sig')
    
    novel_df = final_df[final_df['novelty'] == "NOVEL"]
    novel_df.to_csv(novel_csv, index=False, encoding='utf-8-sig')

    print("\nTURBO-BLAST COMPLETE")
    print("="*70)
    print(f"Audit Trail: {output_csv}")
    print(f"Novel Hits:  {novel_csv} ({len(novel_df)} hits)")
    print("NEXT STEP: Run 'python validation/step3_fasta_prep.py'")
    print("="*70)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PlatyGeno: Turbo-BLAST Validation.")
    parser.add_argument("--input", type=str, default="discovery_hits.csv", help="Discovery CSV to validate")
    parser.add_argument("--output-prefix", type=str, help="Prefix for result files")
    parser.add_argument("--threads", type=int, default=5, help="Number of parallel BLAST threads (default: 5)")
    parser.add_argument("--all", action="store_true", help="Validate ALL sequences (Method A and B)")
    
    args = parser.parse_args()
    
    run_benchmarking_validation(
        input_csv=args.input, 
        output_prefix=args.output_prefix, 
        max_workers=args.threads,
        validate_all=args.all
    )
