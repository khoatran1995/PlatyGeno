import gzip
import shutil
import os
import time

def fast_unzip(accession, folder="data/full_benchmarking"):
    files = [f"{accession}_1.fastq.gz", f"{accession}_2.fastq.gz"]
    
    for filename in files:
        gz_path = os.path.join(folder, filename)
        out_path = os.path.join(folder, filename.replace(".gz", ""))
        
        if not os.path.exists(gz_path):
            print(f"File not found: {gz_path}")
            continue
            
        print(f"Decompressing {filename} -> {os.path.basename(out_path)}...")
        start_time = time.time()
        
        # Using a large 10MB buffer for faster IO on SSDs
        with gzip.open(gz_path, 'rb') as f_in:
            with open(out_path, 'wb', buffering=10*1024*1024) as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        elapsed = time.time() - start_time
        size_gb = os.path.getsize(out_path) / 1024 / 1024 / 1024
        print(f"Finished {filename}. Extracted {size_gb:.1f} GB in {elapsed:.1f}s")

if __name__ == "__main__":
    fast_unzip("SRR23196177")
