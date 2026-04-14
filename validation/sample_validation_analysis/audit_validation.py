import pandas as pd
import numpy as np
from scipy import stats
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
stage1_path = os.path.join(BASE_DIR, "PLG_Stage1_Significance.csv")
stage2_path = os.path.join(BASE_DIR, "PLG_Stage2_Validation.csv")

def thorough_audit():
    print("="*80)
    print("PLATYGENO v1.0.2: DEEP DATA AUDIT & COHESION CHECK")
    print("="*80)

    if not os.path.exists(stage2_path):
        print("[X] Missing Validation Data")
        return

    df2 = pd.read_csv(stage2_path)
    df2['len'] = df2['sequence'].str.len()
    
    # Check for paired features (those having both snippet and assembly)
    p_ids = df2[df2['method'] == 'Precision Snippet']['feature_id'].unique()
    a_ids = df2[df2['method'] == 'Consensus Assembly']['feature_id'].unique()
    paired_ids = set(p_ids).intersection(set(a_ids))
    
    paired_df = df2[df2['feature_id'].isin(paired_ids)].copy()
    
    print(f"\n[PAIRED FEATURE ANALYSIS (N={len(paired_ids)} features)]")
    snippets = paired_df[paired_df['method'] == 'Precision Snippet']
    assemblies = paired_df[paired_df['method'] == 'Consensus Assembly']
    
    print(f"Snippet Avg Length: {snippets['len'].mean():.2f} ± {snippets['len'].std():.2f}")
    print(f"Assembly Avg Length: {assemblies['len'].mean():.2f} ± {assemblies['len'].std():.2f}")
    
    print(f"Snippet Avg Identity: {snippets['blast_identity'].mean():.2f}%")
    print(f"Assembly Avg Identity: {assemblies['blast_identity'].mean():.2f}%")
    
    print(f"Snippet E-v Median: {snippets['e_value'].median():.2e}")
    print(f"Assembly E-v Median: {assemblies['e_value'].median():.2e}")

    # Correlation Analysis
    # Does the report correlation only look at Assemblies? (Since snippets are fixed at 60?)
    # Valid E-values ONLY
    valid_a = assemblies[assemblies['e_value'] > 0].copy()
    valid_a['log_e'] = -np.log10(valid_a['e_value'])
    r, p_val = stats.pearsonr(valid_a['len'], valid_a['log_e'])
    print(f"\n[CORRELATION STUDY (Assemblies Only)]")
    print(f"Pearson r (Length vs -log10E): {r:.4f}")

    # Significance Bins
    def bin_logic(l):
        if l <= 65: return "60 (Snippet Range)"
        if l <= 100: return "80-100"
        return "100+ (Assemblies)"
    
    valid_a['bin'] = valid_a['len'].apply(bin_logic)
    print("\n[SIGNIFICANCE BY LENGTH BIN]")
    print(valid_a.groupby('bin', observed=True)['log_e'].mean())

    # Census verification
    if os.path.exists(stage1_path):
        df1 = pd.read_csv(stage1_path)
        print(f"\n[GLOBAL CENSUS VERIFICATION]")
        print(f"Complexity: {df1['complexity'].mean():.2f} ± {df1['complexity'].std():.2f}")
        print(f"GC%: {df1['gc_content'].mean()*100:.2f}% ± {df1['gc_content'].std()*100:.2f}%")

if __name__ == "__main__":
    thorough_audit()
