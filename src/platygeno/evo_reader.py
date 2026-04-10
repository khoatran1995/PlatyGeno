# Copyright 2026 [Your Name]
# Licensed under the Apache License, Version 2.0 (the "License")

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
from itertools import islice

def read_evo_features(fasta_path, engine, start_idx=0, end_idx=1000, layer=26, score_thres=10.0):
    """
    Function 1: Scans a FASTA/FASTQ range through Evo 2 and maps Layer 26 activations.
    
    Args:
        fasta_path (str): Path to genomic file.
        engine (PlatyGenoEngine): Initialized model engine.
        start_idx (int): Start index for reads.
        end_idx (int): End index for reads.
        layer (int): Evo 2 layer to interpret (default 26).
        score_thres (float): Minimum activation strength to record.
    """
    results = []
    
    # Efficiently skip to start_idx and stop at end_idx
    with open(fasta_path, "r") as handle:
        records = islice(SeqIO.parse(handle, "fasta"), start_idx, end_idx)
        
        print(f"🧬 EvoReader: Analyzing reads {start_idx} to {end_idx}...")
        for i, record in enumerate(tqdm(records, total=(end_idx - start_idx))):
            # Get activations from the SAE
            activations = engine.get_features(str(record.seq))
            
            # Identify indices above the score_thres
            indices = activations.nonzero(as_tuple=True)[0]
            for idx in indices:
                val = activations[idx].item()
                if val >= score_thres:
                    results.append({
                        "read_id": record.id,
                        "feature_id": int(idx.item()),
                        "score": round(val, 4),
                        "sequence": str(record.seq)
                    })
            
    return pd.DataFrame(results)
