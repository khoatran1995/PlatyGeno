import pandas as pd
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
novel_path = os.path.join(BASE_DIR, "PLG_Stage2_Novel_Sequences.csv")
stage1_path = os.path.join(BASE_DIR, "PLG_Stage1_Significance.csv")

def profile_novelty():
    print("="*60)
    print("ANALYSIS 3: NOVEL GENOMIC HITS PROFILE")
    print("="*60)

    if not os.path.exists(novel_path):
        print(f"[ERROR] Missing: {novel_path}")
        return

    df = pd.read_csv(novel_path)
    df['Length'] = df['sequence'].str.len()
    
    # Export Table
    novel_table = df[['feature_id', 'method', 'Length', 'e_value']].copy()
    novel_table.columns = ['Feature ID', 'Assembly Method', 'Length (bp)', 'BLAST E-value']
    
    table_md = novel_table.to_markdown(index=False)
    with open(os.path.join(BASE_DIR, "table_novelty.md"), "w") as f:
        f.write(table_md)

    print(f"\n[SUCCESS] Exported Novelty Table to BASE_DIR")

if __name__ == "__main__":
    profile_novelty()
