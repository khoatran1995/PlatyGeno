# Gene Discovery Validation: Methodology & Audit Trail 🧪

This document provides the technical and scientific context for the **Standard Validation Suite** located in `/validation`. This validation proves the significance of PlatyGeno for unsupervised gene discovery in metagenomic environments.

## 1. Data Sourcing (Provenance)
*   **Dataset**: `sample.fastq`
*   **Origin**: Subsampled from the first 20,000 reads of clinical sample **HSMA33OT_R1** of **The Inflammatory Bowel Disease Multi'omics Database (IBDMDB)**.
*   **Rationale**: The IBDMDB gut metagenome represents a high-impact clinical environment where unsupervised discovery can reveal novel functional sequences associated with direct human health outcomes.

## 2. Processing Engine (The "Golden State")
The validation utilizes the PlatyGeno core engine optimized for detection sensitivity:
*   **Architecture**: Evo 2 7B + Sparse Autoencoders (SAEs).
*   **Scanning Logic**: **Batched Mean-Pooling** across reads to reduce noise and emphasize persistent biological signals.
*   **Significance Filtering**: The **"Padding Filter"** is applied to eliminate artifacts from genomic padding, ensuring that only true biological motifs are surfaced.

## 3. Three-Step Discovery Pipeline
The methodology is executed through three serialized stages found in `/validation`:

### Stage 1: Significance Scanning (`step1_discovery.py`)
*   **Action**: Maps the genomic landscape of the raw dataset using automated SAE interpretation.
*   **Metric**: Activation threshold (8.0+).
*   **Output**: `PLG_Stage1_Significance.csv`. Identifies high-confidence DNA landmarks.

### Stage 2: Database Identification (`step2_blast.py`)
*   **Action**: Performs automated BLAST searches on the candidates identified in Stage 1.
*   **Goal**: Categorize discoveries into "Known Biological Landmarks" (database matches) and "Novel Feature Candidates" (zero-homology).
*   **Output**: `PLG_Stage2_Validation.csv`.

### Stage 3: Discovery Refinement (`step3_fasta_prep.py`)
*   **Action**: Extracts clean FASTA sequences for the top candidates for downstream validation to verify the structure of possible novel, never-before-seen DNA.
*   **Structural Validation**: Prepares sequences for structural modeling (e.g., AlphaFold/ESMFold). A high pLDDT score (>70) in a sequence with no homology match confirms a high-confidence **Novel Discovery**.

## 4. Conclusion
This streamlined methodology moves beyond simple "List Matching." By focusing on a high-complexity clinical gut dataset, PlatyGeno proves its ability to autonomously identify functional genomic features that are clinically relevant, regardless of whether they have been previously cataloged.

---
Developed by **Khoa Tu Tran** for the PlatyGeno Gene Discovery Project.
