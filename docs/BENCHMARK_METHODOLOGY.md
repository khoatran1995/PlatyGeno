# Viral Discovery Benchmark: Methodology & Audit Trail 🧪

This document provides the full technical and scientific context for the **Viral Discovery Showcase** (`examples/viral_discovery_showcase.py`). This benchmark was designed to prove the clinical and biological significance of PlatyGeno for unsupervised gene discovery.

## 1. Data Sourcing (Provenance)
*   **Dataset**: `small_virome.fastq.gz`
*   **Origin**: Subsampled from a validated clinical viral metagenome used in the *Virome Benchmarking Study* (Beard Group).
*   **Rationale**: We selected a clinical virome because viral "Dark Matter" often lacks high-accuracy labels, making it the perfect test case for Evo 2's unsupervised foundation features.

## 2. Processing Pipeline (Three-Phase Discovery)

### Phase 1: Unsupervised Feature Scanning
*   **Engine**: PlatyGeno (Evo 2 7B + Goodfire SAE).
*   **Hook Location**: Transformer Layer 26 (Biological Motif Layer).
*   **Filtering**: Using the "Rare Signal Filter" to isolate features with low frequency but high activation (threshold = 8.0+). This effectively ignores human background DNA and focuses on unique functional viral markers.

### Phase 2: Translation & ORF Recovery
*   **Strategy**: 6-frame translation of DNA discovery snippets.
*   **Heuristic**: "Longest ORF" (Open Reading Frame) selection. This assumes that the functional peak found by the AI corresponds to the coding region of a protein.

### Phase 3: Structural Validation (The Proof)
*   **Inference Engine**: Meta **ESMFold v1**.
*   **Metric**: **pLDDT** (Predicted Local Distance Difference Test).
*   **Significance**:
    *   **pLDDT > 70**: High confidence. The AI-discovered snippet folds into a physically viable protein structure.
    *   **pLDDT > 90**: Extremely high confidence. Equivalent to experimental accuracy.

## 3. How to Interpret the Final Report
The showcase generates `data/benchmark_results_summary.csv`.
*   **`feature_id`**: The "Concept ID" inside the model.
*   **`dna_len`**: Length of the raw DNA snippet extracted.
*   **`plddt`**: The quality score. A high score here is the objective proof that PlatyGeno found a real, functional gene from raw noise.

---
Developed by **Khoa Tu Tran** for the PlatyGeno Gene Discovery Project.
