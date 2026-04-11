# Copyright 2026 Khoa Tu Tran
# Licensed under the Apache License, Version 2.0 (the "License")
#
# gene_discovery_sample.py
#
# Full 3-Phase Gene Discovery Workflow using PlatyGeno.
#
# PHASE 1 — Discovery Scan:
#   Run the Evo 2 model + SAE over all reads in a FASTQ/FASTA file.
#   Output: A feature report CSV mapping read_id -> feature_id -> activation.
#
# PHASE 2 — Rare Signal Filtering:
#   Apply statistical filters to isolate 'Rare but Powerful' features.
#   These are the biological needle-in-a-haystack signals.
#   Output: Top-N candidate feature IDs and their best-scoring read IDs.
#
# PHASE 3 — Precise Gene Snippet Extraction:
#   For each winning read, re-run the deep token-level scan.
#   Isolate the exact 60bp nucleotide window responsible for the high activation.
#   Output: A final results table saved as a CSV.

import os
import torch
import _codecs
import pandas as pd
from Bio import SeqIO

# --- Skeleton Key: torch serialization safety bypass for SAE weights ---
torch.serialization.add_safe_globals([_codecs.encode])

from platygeno.core import PlatyGenoEngine
from platygeno.evo_reader import read_evo_features
from platygeno.mapper import (
    find_rare_needle_signals,
    get_best_reads_for_features,
    extract_precise_gene_code,
)

# =============================================================================
# CONFIGURATION
# =============================================================================
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FASTQ_FILE   = os.path.join(base_dir, "data", "sample.fastq")
PHASE1_CSV   = os.path.join(base_dir, "data", "phase1_feature_report.csv")
RESULTS_CSV  = os.path.join(base_dir, "data", "gene_snippets_top10.csv")

SCAN_LIMIT   = 4000   # Max reads to process in Phase 1
TOP_N        = 10     # Number of candidate genes to extract in Phase 3
WINDOW_SIZE  = 60     # Size (bp) of the extracted gene snippet

# =============================================================================
# SETUP
# =============================================================================
if not os.path.exists(FASTQ_FILE):
    raise FileNotFoundError(
        f"⚠️  FASTQ file not found at: {FASTQ_FILE}\n"
        f"    Please place your input file at the path above."
    )

# Auto-detect format from file extension
ext = os.path.splitext(FASTQ_FILE)[1].lower()
FILE_FORMAT = "fastq" if ext in [".fastq", ".fq"] else "fasta"

# =============================================================================
# PHASE 1: DISCOVERY SCAN
# =============================================================================
print("\n" + "="*60)
print("🚀 PHASE 1: Starting PlatyGeno Engine and Discovery Scan")
print("="*60)

engine = PlatyGenoEngine(model_name='evo2_7b')

print(f"\n📡 Scanning up to {SCAN_LIMIT} reads from: {FASTQ_FILE}")
phase1_report = read_evo_features(FASTQ_FILE, engine, limit=SCAN_LIMIT)

if phase1_report.empty:
    raise RuntimeError("❌ Phase 1 returned no results. Check your input file or model setup.")

phase1_report.to_csv(PHASE1_CSV, index=False)
print(f"\n✅ Phase 1 Complete — {len(phase1_report)} feature activations recorded.")
print(f"   Report saved to: {PHASE1_CSV}")

# =============================================================================
# PHASE 2: RARE SIGNAL FILTERING
# =============================================================================
print("\n" + "="*60)
print("🔬 PHASE 2: Filtering for Rare High-Activation Features")
print("="*60)

# Find the top-N 'Rare Needle' feature candidates
rare_candidates = find_rare_needle_signals(phase1_report, freq_max=40, top_n=TOP_N)

print(f"\n🎯 Top {TOP_N} Rare Feature Candidates:")
print(rare_candidates.to_string(index=False))

# Get the single best (highest activation) read for each candidate feature
winning_feature_ids = rare_candidates['feature_id'].tolist()
winner_reads = get_best_reads_for_features(phase1_report, winning_feature_ids)

print(f"\n🏆 Absolute Winner Reads (one per feature):")
print(winner_reads[['feature_id', 'read_id', 'activation']].to_string(index=False))

# =============================================================================
# PHASE 3: PRECISE GENE SNIPPET EXTRACTION
# =============================================================================
print("\n" + "="*60)
print("🧬 PHASE 3: Deep Token Scan — Extracting Precise Gene Snippets")
print("="*60)

# Build a lookup index: read_id -> DNA sequence from the original file
# We only load the reads we need, not the entire file.
print(f"\n📖 Building sequence index from: {FASTQ_FILE}")
winning_read_ids = set(winner_reads['read_id'].tolist())

# Note: read_id in the report was sanitized (/ | : -> _), but the original
# record.id in the file is raw. We match on the sanitized version.
sequence_index = {}
for record in SeqIO.parse(FASTQ_FILE, FILE_FORMAT):
    safe_id = record.id.replace("/", "_").replace("|", "_").replace(":", "_")
    if safe_id in winning_read_ids:
        sequence_index[safe_id] = str(record.seq)

print(f"   Found {len(sequence_index)} matching sequences in file.")

# Run deep extraction for each winner
final_results = []
for _, row in winner_reads.iterrows():
    feature_id = int(row['feature_id'])
    read_id    = row['read_id']
    phase1_act = row['activation']

    dna_seq = sequence_index.get(read_id)
    if dna_seq is None:
        print(f"   ⚠️  Skipping read '{read_id}': not found in file.")
        continue

    print(f"\n   🔍 Extracting snippet for Feature #{feature_id} from read '{read_id}'...")
    gene_snippet, precise_activation, peak_bp_idx = extract_precise_gene_code(
        engine, dna_seq, feature_id, window_size=WINDOW_SIZE
    )

    final_results.append({
        "rank":               len(final_results) + 1,
        "feature_id":         feature_id,
        "read_id":            read_id,
        "phase1_activation":  round(phase1_act, 4),
        "peak_activation":    round(precise_activation, 4),
        "peak_bp_index":      peak_bp_idx,
        "gene_snippet_60bp":  gene_snippet,
    })
    print(f"   ✅ Gene Snippet: {gene_snippet}")
    print(f"      Peak Activation: {precise_activation:.4f}  |  Peak BP Index: {peak_bp_idx}")

# =============================================================================
# FINAL REPORT
# =============================================================================
results_df = pd.DataFrame(final_results)

print("\n" + "="*60)
print("🎉 DISCOVERY COMPLETE — Final Gene Snippets")
print("="*60)
print(results_df.to_string(index=False))

results_df.to_csv(RESULTS_CSV, index=False)
print(f"\n💾 Results saved to: {RESULTS_CSV}")
