import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
stage2_path = os.path.join(BASE_DIR, "PLG_Stage2_Validation.csv")

def generate_figures():
    print("="*60)
    print("ANALYSIS 4: GENERATING PUBLICATION FIGURES")
    print("="*60)

    if not os.path.exists(stage2_path):
        print(f"[ERROR] Missing: {stage2_path}")
        return

    df = pd.read_csv(stage2_path)
    
    # 1. Figure 1: Significance Gain (E-values)
    # We filter out 0.0 values for log scaling visibility or use rank-based
    plt.figure(figsize=(8, 6))
    
    # Pre-process E-values: Cap at 1e-100 for visualization and handle 0.0
    df['log_eval'] = np.log10(df['e_value'].replace(0, 1e-100))
    
    methods = df['method'].unique()
    data_to_plot = [df[df['method'] == m]['log_eval'] for m in methods]
    
    plt.boxplot(data_to_plot, labels=methods, patch_artist=True, 
                boxprops=dict(facecolor='lightblue', color='blue'),
                medianprops=dict(color='red'))
    
    plt.ylabel("BLAST E-value (log10)")
    plt.title("Statistical Gain: Match Significance by Method")
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    fig1_path = os.path.join(BASE_DIR, "fig_significance.png")
    plt.savefig(fig1_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[SUCCESS] Exported Figure 1: {fig1_path}")

    # 2. Figure 2: Taxonomic Breakdown
    def categorize_hit(title):
        title = str(title).lower()
        if "homo sapiens" in title or "human" in title: return "Host DNA (Human)"
        if any(x in title for x in ["bacteroides", "clostridiales", "faecalibacterium", "vulgatus", "commensal", "gut", "mag:", "sp."]): return "Gut Microbiota"
        if "no hits" in title or "unknown" in title: return "Unclassified"
        return "Other"

    df['category'] = df['top_hit'].apply(categorize_hit)
    tax_counts = df.groupby(['method', 'category']).size().unstack(fill_value=0)
    
    # Plot normalized stacked bar
    tax_norm = tax_counts.div(tax_counts.sum(axis=1), axis=0)
    tax_norm.plot(kind='barh', stacked=True, figsize=(10, 6), colormap='viridis')
    
    plt.title("Taxonomic Relevance: Clinical Sample Context")
    plt.xlabel("Proportion of Discovery Hits")
    plt.legend(title="Biological Origin", bbox_to_anchor=(1.05, 1), loc='upper left')
    
    fig2_path = os.path.join(BASE_DIR, "fig_taxonomy.png")
    plt.savefig(fig2_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[SUCCESS] Exported Figure 2: {fig2_path}")

if __name__ == "__main__":
    generate_figures()
