# Copyright 2026 Khoa Tu Tran
# Licensed under the Apache License, Version 2.0 (the "License")

import os
import gzip
import pandas as pd
from Bio import SeqIO
from itertools import islice
from tqdm import tqdm
import torch

def get_format(file_path):
    """Detects if a file is FASTA or FASTQ, even if compressed."""
    is_gz = file_path.endswith('.gz')
    opener = gzip.open if is_gz else open
    
    try:
        with opener(file_path, 'rt') as f:
            char = f.read(1)
            while char and char.isspace():
                char = f.read(1)
            if char == '>':
                return 'fasta'
            elif char == '@':
                return 'fastq'
    except Exception:
        pass
    return 'fasta' # Default

def read_evo_features(file_path, engine, start=0, stop=4000):
    """
    Reads DNA sequences, extracts SAE features, and returns a report.
    Works for both .fasta and .fastq files, including .gz compressed files.
    Allows range-based scanning (from 'start' to 'stop').
    """
    file_format = get_format(file_path)
    is_gz = file_path.endswith('.gz')
    opener = gzip.open if is_gz else open

    print(f"Detected format: {file_format.upper()} | Data: {os.path.basename(file_path)}")
    print(f"EvoReader: Analyzing reads from index {start} to {stop if stop is not None else 'End'}...")

    discovery_results = []
    
    with opener(file_path, 'rt') as f:
        # Use islice for efficient range-based scanning
        records = islice(SeqIO.parse(f, file_format), start, stop)
        
        for record in tqdm(records, desc="Feature Scanning", unit="read"):
            dna_seq = str(record.seq)
            # Standardize ID for reporting
            safe_id = record.id.replace("/","_").replace("|","_").replace(":","_")
            
            # Extract features (Evo 2 + SAE)
            features = engine.get_read_features(dna_seq)
            
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
