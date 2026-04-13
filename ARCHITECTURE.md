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
*   **Mechanism**: **Mean-Pooling** (Global sequence averaging).
*   **Diversity**: **Zero-Gate Discovery** (Captures ALL active biological signals).
*   **Science**: Filters by `activation` strength (Default: 1.0) and `rarity`.

### 4. The Assembly Layer (Boundary Detection)
*   **Component**: `src/platygeno/workflow.py` & `src/platygeno/mapper.py`.
*   **Function**: Precision extraction of gene boundaries and consensus contig assembly.

---

## 🏗️ Modular & Independent Design

PlatyGeno is built with a strictly decoupled architecture to ensure maximum flexibility for researchers:

### The Independent Core (`src/platygeno`)  
The primary Python package is **completely independent**. It does not require any files from the `validation/` folder to run. You can install it via `pip install -e .` and use `platygeno.discover_genes()` in any script or Jupyter notebook.

### The Validation Suite (`validation/`)  
These are **optional research utilities** designed for Ph.D.-grade verification. They provide "one-touch" automation for BLASTing (via **`step2_blast.py`**) and preparing sequences for AlphaFold, but they are not required for the core engine to function. 
---

## 📦 Version 1.0 Milestone [COMPLETED]

The repository has been successfully upgraded to **Version 1.0**, featuring:
*   **Hardened Core Engine**: Recursive OOM protection and sequence alignment stability.
*   **Unified API**: Simplified `discover_genes()` interface for researchers and notebooks.
*   **Automated Biology**: Integrated Annotation Engine with Found vs. Unknown labeling.
*   **One-Touch Validation**: Ready-to-run 20k benchmark pipeline for gut metagenomes.
*   **Premium Reporting**: Dark-mode HTML dashboards for result visualization.

---

## 🧪 Benchmark Reference: Gut Sample (IBD-MDB)

- **Input**: `data/gut_sample.fastq` (or similar).
- **Goal**: Autonomous discovery of novel features (e.g., Feature 7393).
- **Validation**: BLAST novelty check + AlphaFold structural confidence.
