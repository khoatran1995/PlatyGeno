# PlatyGeno Documentation 🧬

PlatyGeno is a high-performance Python package for **Unsupervised Biological Significance Mapping**. It leverages the Evo 2 foundation model and Sparse Autoencoders to identify genomic landmarks directly from raw sequence data.

---

## 🚀 Quickstart

The most common way to use PlatyGeno is through the automated discovery pipeline:

```python
import platygeno

# Complete Discovery: Scan, Pool, Extract, and Annotate
results = platygeno.discover_genes(
    input_path="data/sample.fastq",
    scan_start=0,
    scan_end=5000,
    min_activation=5.0
)

# View discovered biological features
print(results[['feature_id', 'feature_name', 'activation', 'sequence']])
```

---

## 📚 API Reference (High-Level)

### `platygeno.discover_genes()`
The primary end-to-end interface for genomic discovery missions.

**Parameters:**
- `input_path` (str): Path to FASTQ/FASTA file.
- `engine` (PlatyGenoEngine, optional): Pre-initialized engine instance.
- `scan_start` (int): Starting read index. Default: `0`.
- `scan_end` (int): Ending read index. Set to `None` for full file scan.
- `top_n` (int): Number of features to report. Default: `-1` (All significant hits).
- `top_pct` (float, optional): Return only the top percentage of outliers (e.g., `0.01`).
- `window_size` (int): Genetic snippet extraction window (bp). Default: `60`.
- `min_overlap` (int): Minimum base-pair overlap for landmark assembly. Default: `20`.
- `min_activation` (float): Sensitivity threshold. Default: `1.0`.
- `rel_freq_max` (float): Rarity cap. Set to `1.0` for panoramic discovery.
- `batch_size` (int): GPU batch size (Higher = Faster). Default: `16`.
- `excluded_features` (list[int], optional): List of Feature IDs to ignore.
- `output_path` (str, optional): Destination CSV path.

---

## 🔬 Core Components (Modular Architecture)

For researchers building custom tools, PlatyGeno exposes its underlying discovery modules:

### `platygeno.read_evo_features()`
Generates a raw significance report by passing DNA through the Evo 2 + SAE layers.
- **`file_path`**: Path to DNA sequence data.
- **`engine`**: Initialized `PlatyGenoEngine`.
- **`start / stop`**: Read index range to analyze.
- **`batch_size`**: GPU throughput scaling.

### `platygeno.find_significant_landmarks()`
The statistical core that filters features based on rarity and signal strength.
- **`df`**: Significance report from `read_evo_features`.
- **`rel_freq_max`**: Maximum relative frequency (1.0 = All, 0.001 = Very Rare).
- **`top_n` / `top_pct`**: Number of outlier features to target.
- **`min_activation`**: Minimum SAE signal strength.
- **`excluded_features`**: List of IDs to suppress.

### `platygeno.extract_precise_gene_code()`
Token-level boundary detection for high-resolution mapping.
- **`engine`**: Initialized engine.
- **`dna_seq`**: Raw DNA string.
- **`feature_id`**: The target SAE landmark ID.
- **`window_size`**: Width of the extracted nucleotide snippet (bp).

### `platygeno.assemble_feature_consensus()`
Greedy assembler that merges reads sharing a common biological feature.
- **`sequences`**: List of DNA strings to merge.
- **`min_overlap`**: Minimum overlap (bp) required for a contig merge.

---

## 💡 Technical Note: Why Feature Counts Change

Researchers may notice the number of biological features shifts during the pipeline (e.g., 105 discovered ➔ 101 mapped ➔ 38 validated). This is an intentional feature of the discovery logic:

1.  **Initial Discovery (SAE Activation)**: The engine first finds every unique Feature ID that fires above your threshold.
2.  **Mapping & Quality Control**: A feature is only "Mapped" if it passes the **Extraction & Assembly** phase. If a feature's underlying DNA sequences are too short, noisy, or fail to form a stable consensus overlap, it is filtered out to ensure results are high-confidence.
3.  **Validation (NCBI/BLAST)**: Further drops occur during validation if you filter by "Novelty." Known features (e.g., ribosomal proteins) may be omitted from the final "Novel Candidates" list.

---

## 🏗 Technical Architecture

PlatyGeno layers an interpretability framework on top of genomic foundation models:

1.  **The Signal (Evo 2)**: We extract activations from **Layer 26** of the **Together AI** Evo 2 7B foundation model.
2.  **The Interpretation (SAE)**: We utilize the **Goodfire** Sparse Autoencoders (specifically the `Layer-26-Mixed` expansion) to translate dense AI signals into 32,768 discrete, human-interpretable biological concepts.
3.  **The Result**: Significant signals are extracted, assembled, and labeled via a zero-reference dictionary.

---

## 📜 Citations

If you use PlatyGeno in your research, please cite this repository and the underlying foundation models:

- **PlatyGeno Engine**: Tran, K. T. (2026). *PlatyGeno: Unsupervised Gene Discovery via Evo 2 and Sparse Autoencoders*.
- **Evo 2 Foundation**: Together AI. (2024). *Evo 2: A 7B DNA Foundation Model*.
- **SAE Interpretability**: Goodfire AI. (2024). *Sparse Autoencoders for Genomic Interpretability (Layer 26 Mixed)*.
