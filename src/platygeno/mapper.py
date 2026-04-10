# Copyright 2026 Khoa Tu Tran
# Licensed under the Apache License, Version 2.0 (the "License")

import os
import pandas as pd
import numpy as np

def find_rare_needle_signals(df, freq_max=40, top_n=10):
    """
    Statistical Filter: Groups activations to find 'Rare but Powerful' features.
    This identifies feature IDs that appear in few reads but have high max scores.
    """
    # 1. Calculate population statistics for every feature detected in Phase 1
    feature_stats = df.groupby('feature_id').agg(
        occurrence_count=('read_id', 'count'),
        max_score=('score', 'max')
    ).reset_index()

    # 2. Filter for Rare Candidates (e.g., appearing in < 1% of reads)
    # This ignores common structural motifs and focuses on unique functional genes.
    candidates = feature_stats[
        (feature_stats['occurrence_count'] > 1) & 
        (feature_stats['occurrence_count'] < freq_max)
    ].sort_values(by='max_score', ascending=False)

    return candidates.head(top_n)

def get_best_reads_for_features(df, winning_feature_ids):
    """
    Filters the Phase 1 report to find the specific 'Absolute Winner' reads 
    that triggered the rare features found in find_rare_needle_signals.
    """
    winners = df[df['feature_id'].isin(winning_feature_ids)]
    
    # We only want the single best read for each unique feature ID
    idx = winners.groupby('feature_id')['score'].idxmax()
    return winners.loc[idx].sort_values(by='score', ascending=False)

def extract_precise_gene_code(engine, dna_seq, feature_id, window_size=60):
    """
    Token-Aware Extraction: Re-runs the model only on the candidate sequences
    to isolate the exact nucleotides responsible for the high score.
    """
    # 1. Run the deep token-level scan from core.py
    # token_acts shape: [Seq_Len, 32768]
    token_acts = engine.get_token_features_deep(dna_seq)
    
    # 2. Extract the activation column for the specific feature of interest
    feature_column = token_acts[:, feature_id].cpu().numpy()
    
    # 3. Identify the 'Peak' token (the exact base pair that fired hardest)
    peak_idx = np.argmax(feature_column)
    precise_score = float(feature_column[peak_idx])
    
    # 4. Slice the DNA string around that peak to get the 'Gene Code'
    start = max(0, peak_idx - (window_size // 2))
    end = min(len(dna_seq), start + window_size)
    
    gene_snippet = dna_seq[start:end]
    
    return gene_snippet, precise_score, peak_idx

def annotate_with_biology(df, mapping_file="/workspace/PlatyGeno/data/layer26_features.csv"):
    """Cross-references discovery hits with biological labels (e.g., 'Promoter')."""
    if not os.path.exists(mapping_file):
        print(f"⚠️ Warning: Mapping file {mapping_file} not found.")
        return df
        
    labels_df = pd.read_csv(mapping_file)
    # Merge discovery hits with labels using the feature_id as the key
    annotated = pd.merge(df, labels_df, on='feature_id', how='left')
    return annotated
