import os
import argparse
import json
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

def get_rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def scan_chunk(file_path, start, size, seed_kmers):
    matches = []
    with open(file_path, 'rb') as f:
        f.seek(start)
        if start != 0: f.readline()
        bytes_read = 0
        while bytes_read < size:
            header = f.readline()
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            if not seq: break
            bytes_read += len(header) + len(seq) + len(plus) + len(qual)
            read = seq.strip().decode('ascii').upper()
            if any(km in read for km in seed_kmers):
                matches.append(read)
    return matches

def assemble_strict(seed, matches, min_overlap=50, min_support=2):
    """
    Performs a high-confidence extension.
    Requires 'min_support' different reads to agree on the same extension.
    """
    print(f"\n[PHASE: ASSEMBLE]")
    print(f"Strictness: Overlap={min_overlap}bp, Support={min_support} reads...")
    
    current_seq = seed.upper()
    extended = True
    
    while extended:
        extended = False
        # Extend RIGHT
        suffix = current_seq[-min_overlap:]
        # Get only the portion of the read AFTER the suffix
        candidates = []
        for m in matches:
            idx = m.find(suffix)
            if idx != -1:
                ext_part = m[idx + len(suffix):]
                if ext_part:
                    candidates.append(ext_part)
        
        # Take a 20bp extension for support check
        extensions = [c[:20] for c in candidates if len(c) >= 20]
        
        if extensions:
            counts = Counter(extensions)
            best_ext, support = counts.most_common(1)[0]
            if support >= min_support:
                current_seq += best_ext
                extended = True
                print(f"   -> Extended Right: {len(current_seq)}bp (Support: {support})")
            else:
                print(f"   ! Stopping Right: Support for {best_ext[:5]}... was only {support} (need {min_support})")
        else:
            if candidates:
                print(f"   ! Stopping Right: Candidates found but all were under 20bp.")
            else:
                print(f"   ! Stopping Right: No reads overlap the suffix {suffix[:10]}...")
            
        # Extend LEFT
        prefix = current_seq[:min_overlap]
        candidates = []
        for m in matches:
            idx = m.find(prefix)
            if idx != -1:
                ext_part = m[:idx]
                if ext_part:
                    candidates.append(ext_part)
                    
        extensions = [c[-20:] for c in candidates if len(c) >= 20]
        
        if extensions:
            counts = Counter(extensions)
            best_ext, support = counts.most_common(1)[0]
            if support >= min_support:
                current_seq = best_ext + current_seq
                extended = True
                print(f"   <- Extended Left: {len(current_seq)}bp (Support: {support})")
            else:
                print(f"   ! Stopping Left: Support for ...{best_ext[-5:]} was only {support} (need {min_support})")
        else:
            if candidates:
                print(f"   ! Stopping Left: Candidates found but all were under 20bp.")
            else:
                print(f"   ! Stopping Left: No reads overlap the prefix ...{prefix[-10:]}")
            
    return current_seq

def cmd_scan(args):
    print("="*70)
    print("PLATYGENO DISCOVERY: Phase 1 - High-Confidence Scanning")
    print("="*70)
    
    seed = args.seed.upper()
    seed_rc = get_rev_comp(seed)
    
    # 1. Generate kmers
    k = args.k
    seed_kmers = {seed[i:i+k] for i in range(len(seed) - k + 1)} | \
                 {seed_rc[i:i+k] for i in range(len(seed_rc) - k + 1)}
    
    print(f"Target Seed Length: {len(seed)} bp")
    print(f"Generated {len(seed_kmers)} search kmers (k={k})")
                 
    num_workers = multiprocessing.cpu_count()
    all_hits = []
    fastq_files = [args.r1]
    if args.r2: fastq_files.append(args.r2)
    
    for fastq_file in fastq_files:
        if not os.path.exists(fastq_file):
            print(f"Error: {fastq_file} not found.")
            continue
        print(f"Scanning {fastq_file}...")
        file_size = os.path.getsize(fastq_file)
        chunk_size = file_size // num_workers
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = [executor.submit(scan_chunk, fastq_file, i*chunk_size, chunk_size if i < num_workers-1 else file_size-i*chunk_size, seed_kmers) for i in range(num_workers)]
            for f in futures: all_hits.extend(f.result())
    
    # 2. Extract unique reads
    unique_hits = list(set(all_hits))
    print(f"\nScan Complete: Found {len(all_hits)} matching reads ({len(unique_hits)} unique).")
    
    # 3. Save to JSON
    out_dir = os.path.dirname(args.out)
    if out_dir: os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(unique_hits, f, indent=2)
    
    print(f"Results saved to: {args.out}")
    print("="*70)

def cmd_assemble(args):
    print("="*70)
    print("PLATYGENO DISCOVERY: Phase 2 - Strict Assembly")
    print("="*70)
    
    if not os.path.exists(args.hits):
        print(f"Error: Hits file {args.hits} not found. Run 'scan' first.")
        return
        
    with open(args.hits, "r") as f:
        matches = json.load(f)
        
    # Add reverse complements for assembly
    all_matches = list(set(matches + [get_rev_comp(m) for m in matches]))
    print(f"Loaded {len(matches)} reads. Total effective pool (incl RC): {len(all_matches)}")
    
    # Run extension
    final_seq = assemble_strict(args.seed, all_matches, min_overlap=args.overlap, min_support=args.support)
    
    out_dir = os.path.dirname(args.out)
    if out_dir: os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w") as f:
        f.write(f">feature_7393_strict_fold\n{final_seq}\n")
        
    print(f"\nFinal Assembly: {len(final_seq)} bp")
    print(f"Saved to: {args.out}")
    print("="*70)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PlatyGeno Chimera-Blocker: Two-Step Walker")
    subparsers = parser.add_subparsers(dest="command", required=True)
    
    # Scan Subcommand
    scan_p = subparsers.add_parser("scan", help="Scan FASTQ files for matching reads")
    scan_p.add_argument("--seed", type=str, required=True, help="DNA seed sequence")
    scan_p.add_argument("--r1", type=str, required=True, help="FASTQ R1 file")
    scan_p.add_argument("--r2", type=str, help="FASTQ R2 file (optional)")
    scan_p.add_argument("--out", type=str, default="validation/hits.json", help="Output JSON path")
    scan_p.add_argument("--k", type=int, default=25, help="Kmer length for scanning")
    
    # Assemble Subcommand
    assemble_p = subparsers.add_parser("assemble", help="Assemble hits into a consensus")
    assemble_p.add_argument("--seed", type=str, required=True, help="DNA seed sequence")
    assemble_p.add_argument("--hits", type=str, required=True, help="Input JSON path from scan")
    assemble_p.add_argument("--out", type=str, default="validation/strict_consensus.fasta", help="Output FASTA path")
    assemble_p.add_argument("--overlap", type=int, default=50, help="Minimum overlap for extension")
    assemble_p.add_argument("--support", type=int, default=2, help="Minimum read support for extension")
    
    args = parser.parse_args()
    
    if args.command == "scan":
        cmd_scan(args)
    elif args.command == "assemble":
        cmd_assemble(args)
