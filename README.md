# PlatyGeno 🧬
**Unsupervised Gene Discovery via Evo 2 & Sparse Autoencoders**

[![PyPI version](https://img.shields.io/pypi/v/platygeno.svg)](https://pypi.org/project/platygeno/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

PlatyGeno is a professional Python package designed to interpret the **Evo 2 genomic foundation model**. It bridges the gap between AI interpretability and biological discovery by identifying functional genomic motifs (promoters, enhancers, coding sequences) directly from raw sequence data without requiring labels.

---

## 🔬 Clinical Discovery Benchmark

PlatyGeno includes a professional, PhD-grade benchmarking pipeline designed for clinical metagenome samples (e.g., from the **IBDMDB**). To replicate our high-resolution discovery run (**Top 100 Features**), use the optimized **Cost-Saving Workflow**:

1. **GPU Discovery Phase (RunPod)**:
   ```bash
   python benchmarks/discovery_phase.py
   ```
2. **Local Validation Phase (Local PC)**:
   ```bash
   python benchmarks/validate_local.py
   ```

---

## ⚙️ Hardware Requirements

PlatyGeno is optimized for high-performance genomic research:

> [!IMPORTANT]
> **GPU Mandatory**: An NVIDIA CUDA-enabled GPU (RTX 3090, 4090, A100, or H100) is required for inference.
> **VRAM**: We recommend **24GB VRAM** for the best experience.

*   **Primary Model**: By default, PlatyGeno utilizes the **Evo 2 7B** model. This provides the ideal balance between biological accuracy and memory efficiency for large-scale discovery.

---

## 🚀 Quick Start

The easiest way to use PlatyGeno is via PyPI. Ensure you are on a GPU-enabled instance (A100/H100/RTX 4090).

```bash
# 1. Install PlatyGeno
pip install platygeno

# 2. Install high-performance GPU kernels
pip install flash-attn --no-build-isolation

# 3. Verify Discovery
platygeno --input sample.fastq --threshold 10.0
```

---

## 📦 Usage Modes

### 1. Command Line Interface (CLI)
Perfect for rapid discovery projects without writing any Python code.

```bash
# Basic discovery (uses default threshold 5.0)
platygeno --input data.fastq

# High-precision scan on a specific sequence range
platygeno --input data.fastq --threshold 12.0 --start 0 --end 5000 --output hits.csv
```

### 2. Library Mode (Python API)
Integrate discovery logic into your own bioinformatics pipelines.

```python
import platygeno

# Single-line discovery pipeline
df = platygeno.discover_genes(
    input_path="sample.fastq",
    scan_end=1000,
    min_activation=8.0
)

# Access precision snippets and assembled contigs
print(df[['method', 'activation', 'sequence']].head())
```

---

## 🧩 The Discovery Workflow

PlatyGeno moves from raw data to biological insight in three distinct scientific phases:

```mermaid
graph LR
    A[Raw Data] --> B[Phase 1: Feature Scanning]
    B --> C[Phase 2: Rare Signal Filtering]
    C --> D[Phase 3: Extraction & Assembly]
    D --> E[Biological Discovery]
```

1.  **Phase 1: Feature Scanning**: Streams sequence data through Evo 2 and records activations of Sparse Autoencoder (SAE) nodes.
2.  **Phase 2: Rare Signal Filtering**: Identifies "Rare Needles"—features with low frequency but high activation strength—to filter out common structural noise.
3.  **Phase 3: Extraction & Assembly**: Isolates precise 60bp snippets at activation peaks and merges overlapping reads into longer functional contigs.

---

## 🏗 Package Anatomy

| Module | Purpose | Key Function |
| :--- | :--- | :--- |
| **`workflow.py`** | **Master Pipeline** | `discover_genes()` — High-level entry point. |
| **`core.py`** | **Model Engine** | `PlatyGenoEngine` — Manages Evo 2 & SAE states. |
| **`mapper.py`** | **Bioinformatics** | `assemble_feature_consensus()` — Greedy assembly. |
| **`evo_reader.py`** | **Data Loading** | `read_evo_features()` — Memory-efficient streaming. |

---

## 📊 Understanding Results

Results are saved as a CSV with the following columns:

*   **`method`**: `Best Snippet` (high precision) or `Assembled Contig` (high context).
*   **`feature_id`**: The SAE index (the AI's internal concept of the biological signal).
*   **`activation`**: The strength of the signal. Higher scores indicate stronger feature presence.
*   **`sequence`**: The isolated DNA sequence ready for BLAST or downstream analysis.

---

## 📚 API Reference

### `platygeno.discover_genes()`
| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `input_path` | `str` | *Req* | Path to the `.fasta`/`.fastq` file. |
| `scan_start` | `int` | `0` | First read index to scan. |
| `scan_end` | `int` | `4000` | Last read index to scan. |
| `min_activation` | `float` | `5.0` | Minimum activation strength. |
| `top_n` | `int` | `10` | The number of rare features to target. |
| `output_path` | `str` | `None` | Path to save CSV results. |

---

## 📜 Primary References

If you use **PlatyGeno** in your research, please cite this package along with the foundational works it is built upon:

**1. PlatyGeno (This Package):**
```bibtex
@software{PlatyGeno2026,
  author = {Khoa Tu Tran},
  title = {PlatyGeno: Unsupervised Gene Discovery via Evo 2 & Sparse Autoencoders},
  url = {https://github.com/khoatran1995/PlatyGeno},
  year = {2026}
}
```

**2. Evo 2 (Foundation Model):**
> Arc Institute. (2026). **Genome modeling and design across all domains of life with Evo 2**. *Nature*. [DOI: 10.1038/s41586-026-10176-5]

**3. Evo 2 Interpretability & SAEs:**
> Deng, M., et al. (2025). **Interpreting Evo 2: Arc Institute's Next-Generation Genomic Foundation Model**. *Goodfire Research*. [DOI: 10.5281/zenodo.14895891]

---

## Acknowledgements
Developed by **[Khoa Tu Tran](https://github.com/khoatran1995)**. This project leverages the **[Evo 2](https://github.com/arcinstitute/evo2)** model by the **[Arc Institute](https://arcinstitute.org)** and Sparse Autoencoder architectures by **[Goodfire](https://goodfire.ai)**.

### 📊 Dataset Credits
The clinical benchmarking data used in this project is a subsample of **HSMA33OT_R1** from **The Inflammatory Bowel Disease Multi'omics Database (IBDMDB)**. We thank the IBDMDB investigators for making this high-impact clinical data publicly available for research.
*   **Source**: [ibdmdb.org](https://ibdmdb.org/)
