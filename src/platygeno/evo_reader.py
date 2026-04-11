# Copyright 2026 Khoa Tu Tran
# Licensed under the Apache License, Version 2.0 (the "License")

import os
import pandas as pd
from Bio import SeqIO
from itertools import islice
from tqdm import tqdm
import torch

def read_evo_features(file_path, engine, start=0, stop=4000):
    """
    Reads DNA sequences, extracts SAE features, and returns a report.
    Works for both .fasta and .fastq files.
    Allows range-based scanning (from 'start' to 'stop').
    """
    # 1. Auto-detect file format by peeking at the first non-whitespace character
    file_format = "fasta" # Default
    try:
        with open(file_path, "r") as f:
            char = f.read(1)
            while char and char.isspace():
                char = f.read(1)
            if char == "@":
                file_format = "fastq"
            elif char == ">":
                file_format = "fasta"
            else:
                # Fallback to extension if character is ambiguous
                ext = os.path.splitext(file_path)[1].lower()
                file_format = "fastq" if ext in [".fastq", ".fq"] else "fasta"
    except Exception as e:
        print(f"⚠️ Warning: Format detection peek failed ({e}). Falling back to FASTA.")

    print(f"📡 Detected format: {file_format.upper()} | Data: {os.path.basename(file_path)}")
    
    # 2. Use islice to respect the range [start, stop)
    record_iterator = islice(SeqIO.parse(file_path, file_format), start, stop)
    
    discovery_results = []

    # 3. Processing Loop
    print(f"🧬 EvoReader: Analyzing reads from index {start} to {stop}...")
    total_val = (stop - start) if stop is not None else None
    
    for record in tqdm(record_iterator, total=total_val, desc=f"🧬 EvoReader: Analyzing reads"):
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
