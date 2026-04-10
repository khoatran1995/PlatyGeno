def get_important_reads(df, min_strength=15.0, top_n=10):
    """Function 2: Filter the report for high-impact discovery."""
    return df[df['strength'] >= min_strength].sort_values(by='strength', ascending=False).head(top_n)

def extract_gene_parts(important_df, engine):
    """Function 3: Map activation back to specific sequence segments."""
    # Logic to identify which part of the DNA triggered the feature
    # (Attribution/Saliency logic goes here)
    pass
