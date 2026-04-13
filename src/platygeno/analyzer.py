import os
import pandas as pd

def generate_html_report(results_df, output_path="reports/discovery_report.html"):
    """
    Generates a premium, dark-mode HTML dashboard for genomic discovery results.
    """
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    
    # 1. Prepare Summary Stats
    total_hits = len(results_df)
    unknown_hits = 0
    if 'feature_name' in results_df.columns:
        unknown_hits = len(results_df[results_df['feature_name'] == "Unknown"])
    
    # 2. Build Table Rows
    rows = ""
    for _, row in results_df.iterrows():
        name = row.get('feature_name', 'Unknown')
        role = row.get('biological_role', 'Unknown')
        act = row.get('activation', 0)
        fid = row.get('feature_id', '?')
        length = row.get('length', 0)
        
        row_html = f"""
        <tr>
            <td><span class="badge {'unknown' if name == 'Unknown' else 'found'}">{name if name != 'Unknown' else 'Unknown'}</span></td>
            <td>{fid}</td>
            <td class="role">{role}</td>
            <td class="source">{row.get('source', 'Unknown')}</td>
            <td class="act">{act:.2f}</td>
            <td>{length} bp</td>
        </tr>
        """
        rows += row_html

    # 3. HTML Template
    # (Style block omitted for brevity, updating the badges in the CSS as well)
    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>PlatyGeno | Discovery Dashboard</title>
        <style>
            body {{ font-family: 'Inter', system-ui, sans-serif; background: #0a0c10; color: #e6edf3; padding: 40px; line-height: 1.6; }}
            .container {{ max-width: 1100px; margin: 0 auto; }}
            h1 {{ font-size: 2.5rem; background: linear-gradient(90deg, #7ee787, #aff5b4); -webkit-background-clip: text; -webkit-text-fill-color: transparent; margin-bottom: 5px; }}
            .stats {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 20px; margin: 30px 0; }}
            .card {{ background: #161b22; border: 1px solid #30363d; border-radius: 12px; padding: 20px; text-align: center; }}
            .card h3 {{ color: #8b949e; font-size: 0.9rem; text-transform: uppercase; margin: 0; }}
            .card p {{ font-size: 2.2rem; font-weight: bold; margin: 10px 0 0 0; }}
            table {{ width: 100%; border-collapse: collapse; background: #161b22; border-radius: 12px; overflow: hidden; border: 1px solid #30363d; }}
            th {{ background: #21262d; text-align: left; padding: 15px; color: #8b949e; font-size: 0.8rem; text-transform: uppercase; }}
            td {{ padding: 15px; border-bottom: 1px solid #30363d; font-size: 0.95rem; }}
            .badge {{ padding: 4px 10px; border-radius: 20px; font-size: 0.75rem; font-weight: bold; text-transform: uppercase; }}
            .badge.found {{ background: rgba(35, 134, 54, 0.2); color: #7ee787; border: 1px solid rgba(126, 231, 135, 0.2); }}
            .badge.unknown {{ background: rgba(163, 113, 247, 0.2); color: #d2a8ff; border: 1px solid rgba(210, 168, 255, 0.2); }}
            .act {{ color: #7ee787; font-weight: bold; font-family: monospace; font-size: 1.1rem; }}
            .role {{ color: #8b949e; font-style: italic; }}
            .source {{ color: #58a6ff; font-size: 0.85rem; font-weight: bold; }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>PlatyGeno Discovery Output</h1>
            <p style="color:#8b949e">Unsupervised Significance Mapping via Evo 2 & Sparse Autoencoders</p>
            
            <div class="stats">
                <div class="card"><h3>Total Landmarks</h3><p>{total_hits}</p></div>
                <div class="card"><h3>Unknown</h3><p style="color:#d2a8ff">{unknown_hits}</p></div>
                <div class="card"><h3>Found</h3><p style="color:#7ee787">{total_hits - unknown_hits}</p></div>
            </div>

            <table>
                <thead>
                    <tr>
                        <th>Classification</th>
                        <th>Feature ID</th>
                        <th>Biological Hypothesis</th>
                        <th>Source</th>
                        <th>Significance peak</th>
                        <th>Length</th>
                    </tr>
                </thead>
                <tbody>
                    {rows}
                </tbody>
            </table>
        </div>
    </body>
    </html>
    """
    
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)
    
    print(f"📊 Premium HTML Report generated: {output_path}")

if __name__ == "__main__":
    # Test with dummy data if run directly
    df = pd.DataFrame([{
        'feature_id': 7393, 'feature_name': 'Unknown', 'biological_role': 'Novel Discovery', 'activation': 12.5, 'length': 450
    }])
    generate_html_report(df)
