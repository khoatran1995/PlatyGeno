# Copyright 2026 Khoa Tu Tran
# Licensed under the Apache License, Version 2.0 (the "License")

import os
import pandas as pd
import numpy as np

def find_rare_needle_signals(df, rel_freq_max=0.0005, top_n=10, min_activation=5.0, total_population=None):
    """
    Scale-Aware Statistical Filter: Finds 'Rare but Powerful' features.
    
    Args:
        rel_freq_max: Maximum allowed percentage of reads containing the feature (default 0.05%).
        total_population: Total number of reads in the sample (used to calculate relative rarity).
    """
    # 1. Calculate population statistics
    feature_stats = df.groupby('feature_id').agg(
        occurrence_count=('read_id', 'count'),
        max_score=('activation', 'max')
    ).reset_index()

    # 2. Determine Dynamic Frequency Cap
    # If total_population isn't provided, we estimate it from the unique reads in the report
    total_processed = total_population if total_population else df['read_id'].nunique()
    dynamic_freq_max = max(2, int(total_processed * rel_freq_max))
    
    print(f"   Scale-Aware Filter: Total Population = {total_processed} | Rarity Cap = {dynamic_freq_max} reads ({rel_freq_max*100:.3f}%)")

    # 3. Filter for candidates
    candidates = feature_stats[
        (feature_stats['occurrence_count'] >= 2) & # Min 2 to avoid noise/errors
        (feature_stats['occurrence_count'] <= dynamic_freq_max) &
        (feature_stats['max_score'] >= min_activation)
    ].copy()
    
    # Add Rarity Percentage to metadata
    candidates['rarity_pct'] = (candidates['occurrence_count'] / total_processed) * 100
    
    return candidates.sort_values(by='max_score', ascending=False).head(top_n)

def get_best_reads_for_features(df, winning_feature_ids):
    """
    Filters the Phase 1 report to find the specific 'Absolute Winner' reads 
    that triggered the rare features found in find_rare_needle_signals.
    """
    winners = df[df['feature_id'].isin(winning_feature_ids)]
    
    # We only want the single best read for each unique feature ID
    idx = winners.groupby('feature_id')['activation'].idxmax()
    return winners.loc[idx].sort_values(by='activation', ascending=False)

def extract_precise_gene_code(engine, dna_seq, feature_id, window_size=60):
    """
    Token-Aware Extraction: Re-runs the model only on the candidate sequences
    to isolate the exact nucleotides responsible for the high activation.
    """
    # 1. Run the deep token-level scan from core.py
    # token_acts shape: [Seq_Len, 32768]
    token_acts = engine.get_token_features_deep(dna_seq)
    
    # 2. Extract the activation column for the specific feature of interest
    feature_column = token_acts[:, feature_id].cpu().numpy()
    
    # 3. Identify the 'Peak' token (the exact base pair that fired hardest)
    peak_idx = np.argmax(feature_column)
    precise_activation = float(feature_column[peak_idx])
    
    # 4. Slice the DNA string around that peak to get the 'Gene Code'
    start = max(0, peak_idx - (window_size // 2))
    end = min(len(dna_seq), start + window_size)
    
    gene_snippet = dna_seq[start:end]
    
    return gene_snippet, precise_activation, peak_idx

def _get_overlap(s1, s2, min_len=20):
    """Calculates the suffix-prefix overlap between two strings."""
    max_len = min(len(s1), len(s2))
    for length in range(max_len, min_len - 1, -1):
        if s1.endswith(s2[:length]):
            return length
    return 0

def assemble_feature_consensus(sequences, min_overlap=20):
    """
    Greedy Assembler: Merges DNA sequences sharing a common feature into a contig.
    Starts with the highest activation read and extends it using overlaps.
    """
    if not sequences:
        return ""
    if len(sequences) == 1:
        return sequences[0]

    # Use the first sequence (highest activation) as the seed
    contig = sequences[0]
    unused = list(sequences[1:])
    
    changed = True
    while changed and unused:
        changed = False
        best_overlap = 0
        best_idx = -1
        mode = None # 'suffix' or 'prefix'

        for i, seq in enumerate(unused):
            # 1. Check if seq is a suffix of contig
            ov_s = _get_overlap(contig, seq, min_overlap)
            if ov_s > best_overlap:
                best_overlap = ov_s
                best_idx = i
                mode = 'suffix'
            
            # 2. Check if contig is a suffix of seq (seq is a prefix)
            ov_p = _get_overlap(seq, contig, min_overlap)
            if ov_p > best_overlap:
                best_overlap = ov_p
                best_idx = i
                mode = 'prefix'

        if best_idx != -1:
            if mode == 'suffix':
                contig = contig + unused[best_idx][best_overlap:]
            else:
                contig = unused[best_idx] + contig[best_overlap:]
            unused.pop(best_idx)
            changed = True
            
    return contig

def annotate_with_biology(df, mapping_file=None):
    """Cross-references discovery hits with biological labels (e.g., 'Promoter')."""
    if mapping_file is None:
        # Default to the data directory in the project root
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        mapping_file = os.path.join(base_dir, "data", "layer26_features.csv")
        
    if not os.path.exists(mapping_file):
        print(f"⚠️ Warning: Mapping file {mapping_file} not found.")
        return df
        
    labels_df = pd.read_csv(mapping_file)
    # Merge discovery hits with labels using the feature_id as the key
    annotated = pd.merge(df, labels_df, on='feature_id', how='left')
    return annotated
