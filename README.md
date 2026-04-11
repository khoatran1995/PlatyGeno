# PlatyGeno 🧬
**Unsupervised Gene Discovery via Evo 2 & Sparse Autoencoders**

PlatyGeno is a professional Python package designed to interpret the **Evo 2 genomic foundation model**. It bridges the gap between AI interpretability and biological discovery by identifying functional genomic motifs (promoters, enhancers, coding sequences) without requiring labels.

---

## ⚙️ Hardware Requirements

PlatyGeno is built for high-performance genomic research:

*   **GPU Required**: An NVIDIA CUDA-enabled GPU (e.g., RTX 3090/4090, A100, or H100) is **mandatory** for inference.
*   **Optimal Feasibility**: By default, PlatyGeno utilizes the **Evo 2 7B** model. This provides the ideal balance between biological accuracy and memory efficiency, enabling large-scale discovery on single-GPU nodes.
*   **VRAM**: A minimum of **24GB VRAM** is recommended for stable performance.

---

## 🚀 Quick Start (Production)

The easiest way to run PlatyGeno is via PyPI on a GPU-enabled instance (A100, H100, or RTX 3090/4090).

```bash
# 1. Install from PyPI
pip install platygeno

# 2. Install high-performance GPU kernels
pip install flash-attn --no-build-isolation

# 3. Run the CLI
platygeno --input sample.fastq --threshold 10.0
```

---

## 📦 Usage

### 1. Command Line Interface (CLI)
PlatyGeno comes with a built-in CLI for rapid discovery without writing code.

```bash
platygeno --input your_data.fastq --start 0 --end 1000 --output results.csv
```

### 2. Python API
Integrate PlatyGeno into your own bioinformatics pipelines.

```python
import platygeno

# Discover genes in a single line
df = platygeno.discover_genes(
    input_path="sample.fastq",
    scan_start=0,
    scan_end=1000,
    min_activation=8.5,
    output_path="discovery_report.csv"
)

print(df.head())
```

---

## 🏗 Package Anatomy
PlatyGeno is designed to be modular. You can use the high-level pipeline or build custom workflows using the core modules.

| Module | Purpose | Key Function |
| :--- | :--- | :--- |
| **`workflow.py`** | **Automated Pipeline** | `discover_genes()` — The one-stop-shop for discovery. |
| **`core.py`** | **Model & SAE Engine** | `PlatyGenoEngine` — Manages high-mem model loading. |
| **`mapper.py`** | **Signal Analysis** | `assemble_feature_consensus()` — Greedy contig assembly. |
| **`evo_reader.py`** | **Data Streaming** | `read_evo_features()` — High-speed file scanning. |

### 🛠 Module Breakdown

#### 1. `workflow.py` (The Master Pipeline)
The high-level controller that automates the entire discovery process. It is the easiest way to use the library.
*   **How to use:** 
    ```python
    df = platygeno.discover_genes(
        input_path="sample.fastq", 
        scan_start=0, 
        scan_end=5000, 
        top_n=10, 
        min_activation=5.0, 
        window_size=60, 
        min_overlap=20
    )
    ```

#### 2. `core.py` (The Model Engine)
Manages the connection between the **Evo 2 Foundation Model** and the **Sparse Autoencoder (SAE)**. It handles GPU allocation and layer-26 feature extraction.
*   **How to use:** 
    ```python
    engine = platygeno.PlatyGenoEngine(model_name='evo2_7b', device='cuda')
    ```

#### 3. `mapper.py` (The Bioinformatic Assembly)
Contains the greedy assembly algorithm and signal filtering logic. It "connects" short reads into long, biologically significant contigs based on shared SAE features.
*   **How to use:** 
    ```python
    contig = platygeno.assemble_feature_consensus(overlapping_reads, min_overlap=20)
    ```

#### 4. `evo_reader.py` (The File Streamer)
Provides memory-efficient streaming for massive `.fasta` or `.fastq` files. It uses `islice` to scan specific chunks of a file without loading the whole thing into RAM.
*   **How to use:** 
    ```python
    report = platygeno.read_evo_features("data.fastq", engine, start=0, stop=1000)
    ```

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

## 📜 Citation

If you use **PlatyGeno** in your research, please cite it as follows:

**APA:**
Tran, K. T. (2026). PlatyGeno: Unsupervised Gene Discovery via Evo 2 & Sparse Autoencoders. GitHub. https://github.com/khoatran1995/PlatyGeno

**BibTeX:**
```bibtex
@software{PlatyGeno2026,
  author = {Khoa Tu Tran},
  title = {PlatyGeno: Unsupervised Gene Discovery via Evo 2 & Sparse Autoencoders},
  url = {https://github.com/khoatran1995/PlatyGeno},
  year = {2026}
}
```

---

## Acknowledgements
Developed by **Khoa Tu Tran**. This project leverages the **Evo 2** model by the **Arc Institute** and Sparse Autoencoder architectures by **Goodfire**.
