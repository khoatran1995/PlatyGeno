import pandas as pd
import numpy as np
from scipy import stats
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
stage2_path = os.path.join(BASE_DIR, "PLG_Stage2_Validation.csv")

def compare_techniques():
    print("="*70)
    print("ANALYSIS 2: PRECISION SNIPPET VS. CONSENSUS ASSEMBLY")
    print("="*70)

    if not os.path.exists(stage2_path):
        print(f"[X] Missing: {stage2_path}")
        return

    df = pd.read_csv(stage2_path)
    snippet_df = df[df['method'] == 'Precision Snippet'].copy()
    consensus_df = df[df['method'] == 'Consensus Assembly'].copy()

    snippet_identity = snippet_df['blast_identity'].dropna()
    consensus_identity = consensus_df['blast_identity'].dropna()
    snippet_df['len'] = snippet_df['sequence'].str.len()
    consensus_df['len'] = consensus_df['sequence'].str.len()
    snippet_eval = snippet_df['e_value'].dropna()
    consensus_eval = consensus_df['e_value'].dropna()

    # 1. Statistical Table
    id_u, id_p = stats.mannwhitneyu(snippet_identity, consensus_identity, alternative='two-sided')
    eval_u, eval_p = stats.mannwhitneyu(snippet_eval, consensus_eval, alternative='two-sided')

    stats_data = [
        ["Identity (%)", snippet_identity.mean(), consensus_identity.mean(), id_p],
        ["Sequence Length (bp)", snippet_df['len'].mean(), consensus_df['len'].mean(), "N/A"],
        ["E-value (Median)", snippet_eval.median(), consensus_eval.median(), eval_p]
    ]
    
    stats_df = pd.DataFrame(stats_data, columns=["Metric", "Precision Snippet", "Consensus Assembly", "p-value"])
    stats_md = stats_df.to_markdown(index=False)
    
    with open(os.path.join(BASE_DIR, "table_stats.md"), "w") as f:
        f.write(stats_md)

    # 2. Taxonomic Table
    def categorize_hit(title):
        title = str(title).lower()
        if "homo sapiens" in title or "human" in title: return "Host DNA (Human)"
        if any(x in title for x in ["bacteroides", "clostridiales", "faecalibacterium", "vulgatus", "commensal", "gut", "mag:", "sp."]): return "Gut Microbiota (Target)"
        if "no hits" in title or "unknown" in title: return "Unclassified / Novel"
        return "Other Biological Hit"

    df['category'] = df['top_hit'].apply(categorize_hit)
    tax_summary = df.groupby(['method', 'category']).size().unstack(fill_value=0).reset_index()
    tax_md = tax_summary.to_markdown(index=False)
    
    with open(os.path.join(BASE_DIR, "table_taxonomy.md"), "w") as f:
        f.write(tax_md)

    print(f"\n[SUCCESS] Exported Tables to BASE_DIR")

if __name__ == "__main__":
    compare_techniques()
