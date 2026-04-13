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

    print(f"Detected format: {file_format.upper()} | Data: {os.path.basename(file_path)}")
    print(f"EvoReader: Analyzing reads {start} to {stop if stop is not None else 'End'} (Batch Size: {batch_size})")

    # 1. Initialize data tracking
    all_data = []
    total_processed = 0
    
    with open(file_path, "r") as f:
        # Filter for the requested range
        records = islice(SeqIO.parse(f, file_format), start, stop)
        
        # 2. Process in chunks for high-speed GPU throughput
        file_name = os.path.basename(file_path)
        pbar = tqdm(chunked_iterable(records, batch_size), desc=f"📡 Scanning {file_name}", unit="batch")
        
        for batch in pbar:
            # Filter out sequences that are too short to be biologically significant (or crash SAE)
            valid_batch = [r for r in batch if len(str(r.seq)) >= 32]
            if not valid_batch: continue
            
            dna_strings = [str(r.seq) for r in valid_batch]
            read_ids = [r.id.replace("/","_").replace("|","_").replace(":","_") for r in valid_batch]
            
            # Forward pass through Evo 2 + SAE
            features = engine.get_features(dna_strings)
            
            if features is not None:
                total_processed += len(valid_batch)
                
                # Convert to sparse representation (only keep active neurons)
                # features shape: [batch, 32768]
                vals, inds = torch.topk(features, k=16, dim=-1) # Focus on top 16 signals per read
                
                # Batch processing results for speed
                v_cpu = vals.cpu().numpy()
                i_cpu = inds.cpu().numpy()
                
                for b_idx in range(len(valid_batch)):
                    for f_idx in range(16):
                        act = v_cpu[b_idx, f_idx]
                        if act > 0:
                            all_data.append({
                                "read_id": read_ids[b_idx],
                                "feature_id": int(i_cpu[b_idx, f_idx]),
                                "activation": float(act)
                            })
                            
    # Return both the report and the precise population count for rarity math
    return pd.DataFrame(all_data), total_processed
