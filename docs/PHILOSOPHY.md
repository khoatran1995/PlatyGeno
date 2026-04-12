# Scientific Philosophy: Zero-Reference Significance Mapping

## 1. The Challenge: Homology Bias
Traditional bioinformatics relies on **Homology Search** (e.g., BLAST). This method identifies genes by comparing them to sequences already stored in databases. While effective for known organisms, Homology Search is inherently biased: it cannot identify functional DNA that lacks a known reference.

## 2. The Solution: Deep Genomic Interpretability
PlatyGeno utilizes **Evo 2**, a foundation model trained on objective genomic patterns across all domains of life. Instead of matching DNA to a list, PlatyGeno identifies **Intrinsically Significant** features.

### The Significance Signal
When Evo 2 "reads" a sequence, its internal neurons (Sparse Autoencoders) fire in response to biological concepts. 
- **Activation Strength** is our primary metric. 
- A high activation score indicates that the AI has detected a high-confidence biological feature (e.g., a viral polymerase or a metabolic promoter), regardless of whether that sequence exists in NCBI.

## 3. Reference-Free Discovery Workflow
Under the PlatyGeno "Significance-First" model, the discovery process is inverted:

1.  **Landscape Mapping**: Identifying all significant DNA landmarks in a raw dataset based on AI confidence.
2.  **Structural Validation**: Using AlphaFold to predict the 3D structure of these landmarks. If the DNA matches nothing in the database, but the 3D structure shows a clear fold, a **Novel Discovery** has been made.
3.  **Digital Subtraction**: Optionally filtering for **Rarity** to isolate signals that are non-standard or exclusive to a specific clinical sample.

## 4. Conclusion
PlatyGeno moves bioinformatics away from "List Matching" and toward "Concept Detection." It empowers researchers to find the most important parts of a metagenome purely by observing the "Fluorescence" of biological signals on the AI interpretability manifold.
