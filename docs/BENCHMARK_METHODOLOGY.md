# Clinical Gene Discovery Benchmark: Methodology & Audit Trail 🧪

This document provides the full technical and scientific context for the **Clinical Discovery Showcase** (`examples/viral_discovery_showcase.py`). This benchmark was designed to prove the clinical and biological significance of PlatyGeno for unsupervised gene discovery.

## 1. Data Sourcing (Provenance)
*   **Dataset**: `benchmark_sample.fastq`
*   **Origin**: Subsampled from the clinical sample **HSMA33OT_R1** of **The Inflammatory Bowel Disease Multi'omics Database (IBDMDB)**.
*   **Rationale**: We selected this clinical IBD sample because it represents a high-impact human health environment where unsupervised gene discovery can reveal unique biomarkers and functional genes associated with gut dysbiosis.

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
