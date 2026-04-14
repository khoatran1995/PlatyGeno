import pandas as pd
df = pd.read_csv('PLG_Stage2_Validation.csv')
p = df[df['method'] == 'Precision Snippet'].set_index('feature_id')
a = df[df['method'] == 'Consensus Assembly'].set_index('feature_id')
common = p.index.intersection(a.index)
for fid in common:
    h1 = str(p.loc[fid, 'top_hit'])
    h2 = str(a.loc[fid, 'top_hit'])
    if any(x in h1.lower() for x in ['mag:', 'cand.', 'scaffold', 'fragment', 'unidentified', 'mutant']):
        print(f"{fid} | {h1} | {h2} | {p.loc[fid, 'e_value']:.2e} | {a.loc[fid, 'e_value']:.2e}")
