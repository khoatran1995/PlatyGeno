# Copyright 2026 Khoa Tu Tran
import os
import torch
import _codecs
import pandas as pd
from Bio import SeqIO
from .core import PlatyGenoEngine
from .evo_reader import read_evo_features
from .mapper import (
    find_rare_needle_signals,
    get_best_reads_for_features,
    extract_precise_gene_code,
    assemble_feature_consensus
)

def discover_genes(
    input_path,
    engine=None,
    scan_start=0,
    scan_end=4000,
    top_n=10,
    top_pct=None,
    window_size=60,
    min_overlap=20,
    min_activation=5.0,
    rel_freq_max=0.001,
    output_path=None
):
    """
    End-to-End Discovery Pipeline: Scan -> Filter -> Extract & Assemble.
    
    Args:
        input_path (str): Path to FASTQ/FASTA file.
        engine (PlatyGenoEngine, optional): Existing engine instance.
        scan_start (int): Index of the first read to scan.
        scan_end (int): Index of the last read to scan.
        top_n (int): Number of top rare features to target.
        top_pct (float): Top percentage of candidate features to select.
        window_size (int): Size (bp) of best hit snippets.
        min_overlap (int): Min overlap (bp) for assembly.
        min_activation (float): Min activation to consider a feature as a 'winner'.
        rel_freq_max (float): Maximum allowed percentage of reads containing the feature.
        output_path (str, optional): Save results to this CSV path.
        
    Returns:
        pd.DataFrame: Results containing both Best Snippets and Assembled Contigs.
    """
    # 1. Setup Environment
    torch.serialization.add_safe_globals([_codecs.encode])
    
    if engine is None:
        print("🚀 Initializing PlatyGeno Engine...")
        engine = PlatyGenoEngine()

    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")

    # 2. Phase 1: Feature Scanning
    print(f"📡 Scanning reads {scan_start} to {scan_end if scan_end is not None else 'End'} from {os.path.basename(input_path)}...")
    report, total_scanned = read_evo_features(input_path, engine, start=scan_start, stop=scan_end)
    
    if report.empty:
        print("⚠️ No features detected in initial scan.")
        return pd.DataFrame()

    # 3. Phase 2: Significance Mapping (Zero-Reference)
    # Identify biological landmarks purely from high-activation signals.
    landmarks = find_significant_landmarks(
        report, 
        rel_freq_max=rel_freq_max, 
        top_n=top_n, 
        top_pct=top_pct,
        min_activation=min_activation,
        total_population=total_scanned
    )
    
    if landmarks.empty:
        print(f"⚠️ No biological landmarks met the significance threshold (>={min_activation}).")
        return pd.DataFrame()

    # Create mappings for the final report
    rarity_map = dict(zip(landmarks['feature_id'], landmarks['rarity_pct']))
    count_map = dict(zip(landmarks['feature_id'], landmarks['occurrence_count']))

    # 4. Phase 3: Landmark Extraction & Assembly
    print(f"🧬 Extracting significant snippets and assembling landmarks...")
    ext = os.path.splitext(input_path)[1].lower()
    fmt = "fastq" if ext in [".fastq", ".fq"] else "fasta"
    
    # Pool sequence IDs for the top landmarks
    winning_ids = set()
    for fid in landmarks['feature_id']:
        f_reads = report[report['feature_id'] == fid].sort_values('activation', ascending=False)
        winning_ids.update(f_reads.head(10)['read_id'].tolist())

    # Build sequence index (only for sequences in the winning pool)
    seq_map = {rec.id.replace("/","_").replace("|","_").replace(":","_"): str(rec.seq) 
               for rec in SeqIO.parse(input_path, fmt) if rec.id.replace("/","_").replace("|","_").replace(":","_") in winning_ids}

    final_results = []

    for fid in landmarks['feature_id']:
        f_reads = report[report['feature_id'] == fid].sort_values('activation', ascending=False)
        best_row = f_reads.iloc[0]
        
        # Method A: Significance Point (Precision Extraction)
        rid = best_row['read_id']
        dna = seq_map.get(rid)
        if dna:
            snippet, act, pos = extract_precise_gene_code(engine, dna, fid, window_size=window_size)
            final_results.append({
                "discovery_type": "Significance-Point",
                "method": "Precision Snippet",
                "feature_id": fid,
                "read_id": rid,
                "activation": round(act, 4),
                "occurrence_count": int(count_map.get(fid, 0)),
                "rarity_pct": round(rarity_map.get(fid, 0), 6),
                "length": len(snippet),
                "sequence": snippet
            })

        # Method B: Landmark Assembly (Cross-Read Consensus)
        top_reads = f_reads.head(10)['read_id'].tolist()
        pool = [seq_map[r] for r in top_reads if r in seq_map]
        
        if len(pool) > 1:
            contig = assemble_feature_consensus(pool, min_overlap=min_overlap)
            final_results.append({
                "discovery_type": "Landmark-Assembly",
                "method": "Consensus Assembly",
                "feature_id": fid,
                "read_id": f"Consensus (N={len(pool)})",
                "activation": round(best_row['activation'], 4),
                "occurrence_count": int(count_map.get(fid, 0)),
                "rarity_pct": round(rarity_map.get(fid, 0), 6),
                "length": len(contig),
                "sequence": contig
            })

    results_df = pd.DataFrame(final_results)

    # 5. Output
    if output_path:
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        results_df.to_csv(output_path, index=False)
        print(f"💾 Results saved to {output_path}")

    return results_df
