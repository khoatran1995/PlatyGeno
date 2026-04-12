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

def chunked_iterable(iterable, size):
    it = iter(iterable)
    while True:
        chunk = list(islice(it, size))
        if not chunk:
            break
        yield chunk

def read_evo_features(file_path, engine, start=0, stop=4000, batch_size=16):
    """
    Reads DNA sequences, extracts SAE features in BATCHES, and returns a report.
    This is the high-speed engine for Significance Mapping.
    """
    file_format = get_format(file_path)
    opener = gzip.open if file_path.endswith('.gz') else open

    print(f"Detected format: {file_format.upper()} | Data: {os.path.basename(file_path)}")
    print(f"EvoReader: Analyzing reads {start} to {stop if stop is not None else 'End'} (Batch Size: {batch_size})")

    discovery_results = []
    scanned_count = 0
    
    with opener(file_path, 'rt') as f:
        # 1. Slice the file to the desired range
        records = islice(SeqIO.parse(f, file_format), start, stop)
        
        # 2. Process in chunks for high-speed GPU throughput
        for batch in tqdm(chunked_iterable(records, batch_size), desc="Batch Scanning", unit="batch"):
            batch_seqs = [str(r.seq) for r in batch]
            batch_ids = [r.id.replace("/","_").replace("|","_").replace(":","_") for r in batch]
            
            # Extract features (Evo 2 + SAE) in one GPU call
            # batch_features shape: [Batch_Size, d_hidden]
            batch_features = engine.get_features(batch_seqs)
            
            if batch_features is not None:
                # 3. Process the results for each sequence in the batch
                for i, features in enumerate(batch_features):
                    scanned_count += 1
                    # Find which features 'fired' (activation > 0)
                    active_indices = torch.where(features > 0)[0]
                    active_values = features[active_indices]
                    
                    for idx, val in zip(active_indices, active_values):
                        discovery_results.append({
                            "read_id": batch_ids[i],
                            "feature_id": int(idx.item()),
                            "activation": float(val.item())
                        })
            else:
                scanned_count += len(batch)

    # Convert to DataFrame
    return pd.DataFrame(discovery_results), scanned_count

    # Convert to DataFrame
    return pd.DataFrame(discovery_results), scanned_count
