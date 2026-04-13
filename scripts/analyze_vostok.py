import os
import torch
import pandas as pd
from tqdm import tqdm
from src.platygeno.core import PlatyGenoEngine

def run_discovery(fastq_path, output_csv, batch_size=32, limit=200000):
    """
    Scans a FASTQ file for 'Novelty Spikes' using Evo2 and SAE activations.
    """
    print(f"🧬 Initializing PlatyGeno Discovery Engine...")
    engine = PlatyGenoEngine()
    
    hits = []
    
    print(f"🔍 Scanning {limit} reads from {fastq_path}...")
    
    with open(fastq_path, "r") as f:
        read_batch = []
        read_ids = []
        count = 0
        
        # tqdm for progress tracking
        pbar = tqdm(total=limit)
        
        while count < limit:
            line = f.readline()
            if not line: break
            
            if line.startswith("@"):
                read_id = line.strip()[1:]
                seq = f.readline().strip()
                f.readline() # +
                f.readline() # qual
                
                read_batch.append(seq)
                read_ids.append(read_id)
                count += 1
                pbar.update(1)
                
                if len(read_batch) >= batch_size:
                    # Process batch
                    features = engine.get_features(read_batch)
                    
                    if features is not None:
                        # Find max activation in each read
                        # features shape: [batch, d_hidden]
                        max_vals, max_indices = torch.max(features, dim=1)
                        
                        for i in range(len(read_batch)):
                            if max_vals[i] > 3.0: # Baseline novelty threshold
                                hits.append({
                                    "read_id": read_ids[i],
                                    "activation": max_vals[i].item(),
                                    "feature_id": max_indices[i].item(),
                                    "sequence": read_batch[i]
                                })
                    
                    read_batch = []
                    read_ids = []

    pbar.close()
    
    if hits:
        df = pd.DataFrame(hits)
        # Sort by activation (highest novelty first)
        df = df.sort_values(by="activation", ascending=False)
        df.to_csv(output_csv, index=False)
        print(f"✅ Discovery Complete! Found {len(hits)} novelty spikes.")
        print(f"Top activation: {df['activation'].max():.2f}")
        print(f"Results saved to: {output_csv}")
    else:
        print("❌ No high-activation features found (Threshold > 3.0).")

if __name__ == "__main__":
    input_fastq = "data/SRR5462529/vostok_200k.fastq"
    output_res = "validation/vostok_discovery.csv"
    
    os.makedirs("validation", exist_ok=True)
    run_discovery(input_fastq, output_res)
