import pandas as pd
import os

# Set working directory to the script's location for robust pathing
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
stage1_path = os.path.join(BASE_DIR, "PLG_Stage1_Significance.csv")
stage2_path = os.path.join(BASE_DIR, "PLG_Stage2_Validation.csv")

def summarize_discovery():
    print("="*60)
    print("ANALYSIS 1: DISCOVERY & VALIDATION SUMMARY")
    print("="*60)

    results = []

    # 1. Stage 1: Significance Summary
    if os.path.exists(stage1_path):
        df1 = pd.read_csv(stage1_path)
        unique_features = df1['feature_id'].nunique()
        metrics = ['activation', 'gc_content', 'complexity', 'length']
        means = df1[metrics].mean()
        
        results.append(["Stage 1: Discovery", len(df1), unique_features, means['activation'], means['gc_content'], means['complexity'], means['length']])
        
        print(f"\n[Stage 1] Unique Features: {unique_features}")
    else:
        print(f"[X] Missing: {stage1_path}")

    # 2. Stage 2: Validation Summary
    if os.path.exists(stage2_path):
        df2 = pd.read_csv(stage2_path)
        unique_validated = df2['feature_id'].nunique()
        results.append(["Stage 2: Validation", len(df2), unique_validated, "N/A", "N/A", "N/A", "N/A"])
        print(f"[Stage 2] Unique Features: {unique_validated}")
    
    # Export to Markdown Table
    if results:
        headers = ["Stage", "Total Units", "Unique Features", "Avg Activation", "Avg GC", "Avg Complexity", "Avg Length"]
        table_df = pd.DataFrame(results, columns=headers)
        table_md = table_df.to_markdown(index=False)
        
        output_file = os.path.join(BASE_DIR, "table_discovery.md")
        with open(output_file, "w") as f:
            f.write(table_md)
        print(f"\n[SUCCESS] Exported Table to: {output_file}")

if __name__ == "__main__":
    summarize_discovery()
