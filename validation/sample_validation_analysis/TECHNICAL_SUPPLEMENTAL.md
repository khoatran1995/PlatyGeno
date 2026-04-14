# Supplemental Methodology: Taxonomic Classification & Data Audit

This document outlines the methodological framework used to categorize the 203 biological features identified by PlatyGeno v1.0.2 in the clinical gut metagenomic dataset.

## 1. Classification Methodology

To facilitate biological interpretation, the "Top Hit" results from the BLAST search in `PLG_Stage2_Validation.csv` were mapped into four distinct functional categories:

### A. Gut Microbiota (Target Habitat)
**Classification Rule**: Identification of key bacterial genera and markers associated with the human intestinal microbiome.
*   **Key Search Terms**: `Bacteroides`, `Phocaeicola`, `Clostridiales`, `Alistipes`, `Lachnospiraceae`, `sp.`, `Gut`, and `MAG:` (Metagenome-Assembled Genomes).
*   **Rationale**: These represent the expected biological signal of the sample habitat.

### B. Host DNA (Human Contamination)
**Classification Rule**: Identification of sequences belonging to the human genome.
*   **Key Search Terms**: `Homo sapiens`, `Human`, `Chromosome`, `Scaffold [with human ID prefixes]`.
*   **Rationale**: These represent "off-target" discoveries identifying host genetic material.

### C. Other Biological Hit (Off-Target Discoveries)
**Classification Rule**: Identification of valid biological hits that do not originate from the gut microbiota or the human host.
*   **Key Search Terms**: Valid BLAST hits ($E \le 10^{-5}$) not matching categories A or B (e.g., environmental contaminants, food-derived DNA, or specific viruses).

### D. Unclassified / Novel
**Classification Rule**: Features that failed to return a high-confidence match in the public BLAST databases.
*   **Threshold**: E-value = **10.0** (or empty/null results).
*   **Rationale**: These candidates (like **Feature 7393**) are flagged for structural validation via AlphaFold2 as potential novel genomic landmarks.

---

## 2. Statistical Analysis Methodology

### 2.1 Pearson Correlation
The Pearson Correlation ($r$) was calculated between **Sequence Length (bp)** and the **Log-Transformed Match Significance** ($-\log_{10}(E)$). 
*   **Population**: Paired collection of 98 biological landmarks ($N=196$ total observations).
*   **Transformation**: Null hits (E=10.0) were excluded from the correlation calculation to focus on the performance gain for identifiable features.

### 2.2 Binning Strategy
Significance bins were calculated as the mean $-\log_{10}(E)$ for features within specific length ranges:
*   **Bin 1 (60bp)**: Isolated Precision Snippets.
*   **Bin 2 (80-100bp)**: Intermediate Consensus Assemblies.
*   **Bin 3 (100bp+)**: Full-Scale Discovery Landmarks.

---

## 3. Computational Audit Trail
*   **Ground Truth Data**: `PLG_Stage2_Validation.csv`
*   **Audit Script**: `audit_validation.py`
*   **Verification Date**: April 15, 2026
