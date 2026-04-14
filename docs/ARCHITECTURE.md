# PlatyGeno Technical Architecture: Scientific Foundation & Engine Logic

This document provides a cohesive deep-dive into the underlying methodology, structural stability, and validation framework of the PlatyGeno discovery engine.

---

## I. Scientific Foundations: Solving the Genomic "Needle-in-a-Haystack" Problem

Traditional genomic analysis often relies on hypothesis-driven searching (e.g., BLAST), which requires researchers to know what they are looking for before they begin. In complex metagenomes containing millions of fragmented reads, this approach creates a discovery bottleneck: novel or divergent biological signals are frequently overlooked because they lack known homology.

PlatyGeno implements a paradigm shift to **Significance-First Scanning**:
*   **Intrinsically Significant Detection**: Instead of matching against a database, the engine identifies biological landmarks based on their "semantic intensity" within the DNA grammar, identified by the Evo 2 foundation model.
*   **Bypassing Reference Bias**: By evaluating DNA based on intrinsic functional signals, PlatyGeno can identify highly significant discoveries (e.g., novel viruses or enzymes) that have never been cataloged in current libraries.

---

## II. The Core Discovery Engine: Layer 26 SAE-Guided Interpretation

PlatyGeno derives its interpretability by layering Sparse Autoencoders (SAEs) on top of the **Evo 2 7B** genomic foundation model.

1.  **Semantic Capture (Layer 26)**: As DNA sequences pass through the model, PlatyGeno intercepts the signal at **Layer 26**. This layer serves as the "functional threshold" where raw base pairs are fully translated into complex biological logic.
2.  **SAE Pinpointing**: We project this dense math through a **32,768-feature Sparse Autoencoder**. This identifies discrete biological "concepts" (e.g., *Sigma-70 promoters*, *metalloproteinases*, *viral structural proteins*) that fire within the sequence context.

---

## III. Denoising & Stability: The Batched Padding Filter

To ensure high-precision discovery in noisy clinical samples, we utilize **Batched Mean-Pooling** (The Padding Filter):

*   **Noise Dilation**: By processing sequences in batches, the engine uses sequence padding to "dilate" and neutralize weaker semantic noise.
*   **Signal Stabilization**: Only high-confidence, persistent biological signals survive the pooling phase, ensuring that the final discovery results represent statistically robust genomic landmarks rather than stochastic artifacts.

---

## IV. Discovery Modalities: Pinpointing Snippets and Reconstructing Consensus

PlatyGeno preserves the discovery signal through two distinct, synergistic modalities:

### 🎯 Precision Snippets: High-Interest DNA Clips
Unlike traditional tools that analyze entire reads, PlatyGeno uses the Layer 26 SAE to pinpoint the exact core of a sequence. The engine identifies the **Peak Activation Token**—the single nucleotide coordinate where a biological concept fires most intensely—and extracts a narrow **clip** centered on that hotspot. This strips away biological "noise," providing a precise beacon of significance.

### 🧬 Consensus Assembly: Complete Sequence Reconstruction
To maximize discovery confidence, PlatyGeno merges multiple Precision Snippets that share a common SAE Feature ID. By aligning these narrow clips like a puzzle, the engine reconstructs **complete genomic sequences** ($L \approx 100bp+$).

> [!IMPORTANT]
> **Performance Highlight**: While both modes are preserved in the results, statistical validation confirms that **Consensus Assembly** provides significantly superior biological significance (lower E-values) and cleaner taxonomic assignments compared to individual snippets.

---

## V. The Validation Audit Trail: From Raw Signals to Clinical Proof

We utilize a 3-stage validation workflow to verify discoveries in complex datasets like the **IBD Metagenomic Database (IBD-MDB)**:

1.  **Stage 1: Significance Scanning**: Maps the genomic landscape using automated SAE analysis to identify high-activation (8.0+) landmarks.
2.  **Stage 2: Database Identification (BLAST)**: Categorizes discoveries. Results typically show ~72% alignment with **Gut Microbiota** (Target) and ~22% with **Host DNA** (Human Context), confirming the tool is locked onto the correct biological sample.
3.  **Stage 3: Structural Discovery (AlphaFold)**: Targets the "Novel" population identified in Stage 2. A high confidence score (pLDDT > 70) in a sequence with no homology confirms a **Novel High-Confidence Discovery** (e.g., Feature 7393).

---
*Developed by **Khoa Tu Tran** for the PlatyGeno Research Project.*
