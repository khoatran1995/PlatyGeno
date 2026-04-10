# Copyright 2026 Khoa Tu Tran
# Licensed under the Apache License, Version 2.0 (the "License")

def get_high_score_reads(df, score_thres=15.0, top_n=None):
    """Function 2: Filter the report for high-impact genomic signals."""
    filtered_df = df[df['score'] >= score_thres]
    sorted_df = filtered_df.sort_values(by='score', ascending=False)
    
    if top_n:
        return sorted_df.head(top_n)
    return sorted_df

def extract_motif_segments(target_df, window_size=50):
    """
    Function 3: Isolate the specific DNA part linked to the feature.
    (This is an 'easy-to-modify' placeholder for future saliency logic)
    """
    # For now, this identifies the read; in the future, we can add
    # sliding-window logic here to pinpoint exactly where the 'score' is highest.
    return target_df[['read_id', 'feature_id', 'score', 'sequence']]
