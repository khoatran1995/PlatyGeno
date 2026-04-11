# PlatyGeno 🧬
**Unsupervised Gene Discovery via Evo 2 & Sparse Autoencoders**

PlatyGeno is a professional Python package designed to interpret the **Evo 2 genomic foundation model**. It bridges the gap between AI interpretability and biological discovery by identifying functional genomic motifs (promoters, enhancers, coding sequences) without requiring labels.

---

## 🚀 Quick Start (RunPod)

The fastest way to run PlatyGeno is on a GPU-enabled instance (A100/H100).

```bash
# 1. Clone & Install
git clone https://github.com/khoatran1995/PlatyGeno.git
cd PlatyGeno
pip install flash-attn --no-build-isolation
pip install -e .

# 2. Run the Automated Pipeline
python examples/all_in_one_discovery.py
```

---

## 🏗 Package Anatomy
PlatyGeno is designed to be modular. You can use the high-level pipeline or build custom workflows using the core modules.

| Module | Purpose | Key Function |
| :--- | :--- | :--- |
| **`workflow.py`** | **Automated Pipeline** | `discover_genes()` — The one-line discovery API. |
| **`core.py`** | **Model & SAE Engine** | `PlatyGenoEngine` — Manages Evo 2 and SAE weights. |
| **`mapper.py`** | **Signal Analysis** | `assemble_feature_consensus()` — Greedy contig builder. |
| **`evo_reader.py`** | **Data Streaming** | `read_evo_features()` — Efficiently scans large FASTQ files. |

---

## 🧬 The Discovery Workflow (3-Phases)

PlatyGeno operates in three distinct scientific phases to move from raw data to biological insight:

### Phase 1: Feature Scanning
The `EvoReader` streams your sequence data through the model. It records when specific "Biological Features" (SAE nodes) are activated.
*   *Scale:* Scans thousands of reads per minute.

### Phase 2: Rare Signal Filtering
The `Mapper` applies statistical filters to find "Rare Needles in the Haystack"—features that don't appear in every read but show powerful activation when they do. This filters out common structural motifs to find unique functional genes.

### Phase 3: Extraction & Assembly
Dual-track processing isolates the exact nucleotide sequence:
1.  **Best Snippet**: A precise 60bp window focused on the activation peak.
2.  **Assembled Contig**: Merges multiple overlapping reads for the same feature to reconstruct longer genomic context (often 100bp+).

---

## 📁 Two Ways to Discover

### 1. All-in-One (Production Mode)
Perfect for standard genomic discovery projects.
```bash
# Script: examples/all_in_one_discovery.py
python examples/all_in_one_discovery.py --threshold 10.0 --start 0 --end 10000
```

### 2. Step-by-Step (Research Mode)
Perfect for building custom workflows or integrating with other tools (like BLAST).
```bash
# Script: examples/step_by_step_discovery.py
python examples/step_by_step_discovery.py
```

---

## 📊 Understanding Results
Results are saved to `data/gene_discovery_results.csv`.

*   **`method`**: `Best Snippet` (high precision) or `Assembled Contig` (high context).
*   **`feature_id`**: The SAE index (the "AI's concept" of the biological signal).
*   **`activation`**: The strength of the signal. Higher is better.
*   **`sequence`**: The DNA sequence ready for BLAST search.

---

## 📚 API Reference

### `platygeno.discover_genes()`
The high-level discovery pipeline used in `all_in_one_discovery.py`.

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `input_path` | `str` | Required | Path to the genomic file (`.fasta`/`.fastq`). |
| `scan_start` | `int` | `0` | First read index to scan. |
| `scan_end` | `int` | `4000` | Last read index to scan. |
| `min_activation` | `float` | `5.0` | Minimum strength of a feature to be considered. |
| `top_n` | `int` | `10` | The number of rare features to target. |
| `window_size` | `int` | `60` | Snippet window size in base-pairs. |
| `min_overlap` | `int` | `20` | Min overlap required for assembly. |
| `output_path` | `str` | `None` | Optional path to save CSV results. |

### `platygeno.PlatyGenoEngine()`
Initializes the Evo 2 and SAE models.

*   `model_name`: Default `'evo2_7b'`.
*   `device`: Auto-detects `'cuda'`.

### `platygeno.assemble_feature_consensus()`
The logic behind the "Assembled Contig" method.
*   **Input**: A list of overlapping DNA strings.
*   **Output**: A single merged consensus string.

---

## Acknowledgements
Developed by **Khoa Tu Tran**. This project leverages the **Evo 2** model by the **Arc Institute** and Sparse Autoencoder architectures by **Goodfire**.
