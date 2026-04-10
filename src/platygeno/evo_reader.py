# Copyright 2026 Khoa Tu Tran
# Licensed under the Apache License, Version 2.0 (the "License")

import os
import pandas as pd
from Bio import SeqIO
from itertools import islice
from tqdm import tqdm
import torch

def read_evo_features(file_path, engine, limit=4000):
    """
    Reads DNA sequences, extracts SAE features, and returns a report.
    Works for both .fasta and .fastq files.
    """
    # 1. Auto-detect file format
    ext = os.path.splitext(file_path)[1].lower()
    file_format = "fastq" if ext in [".fastq", ".fq"] else "fasta"
    
    print(f"📖 Reading {file_format.upper()} file: {file_path}")
    
    # 2. Use islice to respect the limit without loading the whole file
    record_iterator = islice(SeqIO.parse(file_path, file_format), 0, limit)
    
    discovery_results = []

    # 3. Processing Loop
    print(f"🧬 EvoReader: Analyzing up to {limit} reads...")
    for record in tqdm(record_iterator, total=limit):
        safe_id = record.id.replace("/", "_").replace("|", "_").replace(":", "_")
        seq = str(record.seq)
        
        # Get features from our Engine (which uses the SAE)
        features = engine.get_features(seq)
        
        if features is not None:
            # Find which features 'fired' (activation > 0)
            active_indices = torch.where(features > 0)[0]
            active_values = features[active_indices]
            
            for idx, val in zip(active_indices, active_values):
                discovery_results.append({
                    "read_id": safe_id,
                    "feature_id": int(idx.item()),
                    "activation": float(val.item())
                })

    # 4. Convert to DataFrame for the final report
    return pd.DataFrame(discovery_results)
