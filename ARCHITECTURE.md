# PlatyGeno Repository Architecture & Revision Guide

This document serves as the high-level roadmap for the PlatyGeno (EVO2) repository. It outlines the core components, discovery workflows, and the plan for repository-wide revisions.

---

## 🏗️ Technical Architecture

### 1. The Foundation Layer (Evo 2)
*   **Model**: 7B parameter DNA/Protein transformer.
*   **Intelligence**: Understands sequence "grammar" across all domains of life.
*   **Role**: Activation extraction source.

### 2. The Interpretability Layer (SAE Generator)
*   **Component**: `src/platygeno/core.py` -> `SparseAutoencoder`.
*   **Mechanism**: Maps activations to **32,768 discrete concept nodes**.
*   **Role**: Translates math into biological "Features".

### 3. The Discovery Engine (Significance Mapping)
*   **Component**: `src/platygeno/mapper.py` & `src/platygeno/evo_reader.py`.
*   **Function**: Batched scanning and peak detection.
*   **Science**: Filters by `activation` strength and `rarity`.

### 4. The Assembly Layer (Boundary Detection)
*   **Component**: `src/platygeno/workflow.py` & `src/platygeno/mapper.py`.
*   **Function**: Precision extraction of gene boundaries and consensus contig assembly.

---

## 🏗️ Modular & Independent Design

PlatyGeno is built with a strictly decoupled architecture to ensure maximum flexibility for researchers:

### The Independent Core (`src/platygeno`)  
The primary Python package is **completely independent**. It does not require any files from the `validation/` folder to run. You can install it via `pip install -e .` and use `platygeno.discover_genes()` in any script or Jupyter notebook.

### The Validation Suite (`validation/`)  
These are **optional research utilities** designed for Ph.D.-grade verification. They provide "one-touch" automation for BLASTing and preparing sequences for AlphaFold, but they are not required for the core engine to function. 

## 📦 Revision Roadmap (The "Big Picture" Plan)

We will follow this sequence to refine and update the repository to the next version:

### Phase 1: Core Package Hardening (`src/`)
- [ ] **Standardize Batching**: Ensure `core.py` and `workflow.py` handle larger datasets without OOM.
- [ ] **Feature Documentation**: Add internal docstrings for SAE feature IDs.
- [ ] **Error Handling**: Add robust checks for CUDA status and model weights.

### Phase 2: CLI & API Enhancements
- [ ] **Interactive Progress**: Improve the `tqdm` or progress bars in the CLI.
- [ ] **Export Options**: Add direct CSV/FASTA export controls to the CLI.
- [ ] **Unified API**: Simplify `discover_genes()` for easier integration into notebooks.

### Phase 3: Validation Suite Cleanup (`validation/`)
- [ ] **Script Modernization**: Sync `step1-4` scripts with latest package changes.
- [ ] **Sample Automation**: Create a one-touch validation command for the "Gut Sample".

### Phase 4: Documentation & Reporting
- [ ] **Dynamic Reports**: Generate HTML reports with activation peak charts.
- [ ] **Installation Guide**: Clarify GPU requirements and `flash-attn` installation steps.

---

## 🧪 Benchmark Reference: Gut Sample (IBD-MDB)

- **Input**: `data/gut_sample.fastq` (or similar).
- **Goal**: Autonomous discovery of novel features (e.g., Feature 7393).
- **Validation**: BLAST novelty check + AlphaFold structural confidence.
