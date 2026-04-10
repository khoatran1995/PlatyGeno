# Copyright 2026 Khoa Tu Tran
# Licensed under the Apache License, Version 2.0 (the "License")

# Copyright 2026 Khoa Tu Tran
# Licensed under the Apache License, Version 2.0 (the "License")

import os
import pandas as pd

def find_rare_needle_signals(df, freq_max=40, top_n=10):
    """
    Notebook Logic (Cell 29): 
    Groups all 256,000+ activations to find 'Rare but Powerful' features.
    """
    # 1. Group by Feature ID to see how often each biological concept appears
    feature_stats = df.groupby('feature_id').agg(
        occurrence_count=('read_id', 'count'),
        max_score=('score', 'max')
    ).reset_index()

    # 2. Filter for Rare Candidates (Appearing in < freq_max reads)
    # This excludes 'common' DNA patterns and focuses on unique genes.
    candidates = feature_stats[
        (feature_stats['occurrence_count'] > 1) & 
        (feature_stats['occurrence_count'] < freq_max)
    ].sort_values(by='max_score', ascending=False)

    return candidates.head(top_n)

def get_winning_gene_code(df, winning_feature_ids):
    """
    Notebook Logic (Cell 30/32):
    For the top rare features, find the specific reads and DNA strings.
    """
    # Filter the main report for only our 'Needle' features
    winners = df[df['feature_id'].isin(winning_feature_ids)]
    
    # For every unique feature, we only want the single read that had the HIGHEST score
    # This is your 'Absolute Winner DNA'
    idx = winners.groupby('feature_id')['score'].idxmax()
    return winners.loc[idx].sort_values(by='score', ascending=False)

def annotate_with_biology(df, mapping_file="/workspace/PlatyGeno/data/layer26_features.csv"):
    """Cross-references feature IDs with the functional labels from your CSV."""
    if not os.path.exists(mapping_file):
        print(f"⚠️ Warning: {mapping_file} not found. Returning unannotated genes.")
        return df
        
    labels_df = pd.read_csv(mapping_file)
    # Merge discovery hits with biological labels (e.g., 'Promoter', 'Zinc Finger')
    annotated = pd.merge(df, labels_df, on='feature_id', how='left')
    return annotated
