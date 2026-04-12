# PlatyGeno 🧬
**Unsupervised Biological Significance Mapping via Evo 2 & Sparse Autoencoders**

[![PyPI version](https://img.shields.io/pypi/v/platygeno.svg)](https://pypi.org/project/platygeno/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

PlatyGeno is a professional Python package for identifying **genomic landmarks** directly from raw sequence data. By leveraging the **Evo 2 foundation model**, it identifies biologically significant DNA structures (promoters, coding sequences, precise motifs) based purely on AI confidence—**without requiring labels, databases, or BLAST.**

---

![PlatyGeno Banner](file:///C:/Users/Admin/.gemini/antigravity/brain/4f3e50e2-8a73-40d5-900c-be13c257ea2e/platygeno_banner_1775984542734.png)

---

## 🔭 Scientific Philosophy: Zero-Reference Significance

Most bioinformatics tools are designed to find "matches" to known lists. PlatyGeno is a **Reference-Free Microscope** that detects the "Signal" of life itself:

*   **Signals over Samples**: Traditional tools (like BLAST) find genes by comparing them to a library. PlatyGeno "reads" the DNA grammar and detects significance peaks directly. If a sequence is important, the AI will find it—even if it's never been sequenced before.
*   **Significance First**: We prioritize **Activation Strength** (the intensity of the AI's internal response). A high activation score is a "Biological Beacon" that points to a functional region.
*   **Optional Novelty Mining**: Once significant landmarks are identified, researchers can optionally filter for **Rarity** to isolate "Genomic Dark Matter" (novel viruses, exotic enzymes, or extremophiles).

---

## ⚙️ Installation

PlatyGeno requires a CUDA-enabled GPU (RTX 3090, 4090, A100, or H100).

```bash
# 1. Install the core package
pip install platygeno

# 2. Install high-performance GPU kernels (Mandatory for speed)
pip install flash-attn --no-build-isolation

# 3. Development/Editable Install
git clone https://github.com/khoatran1995/PlatyGeno.git
cd PlatyGeno
pip install -e .
```

---

## 🏗️ Simplified Architecture

PlatyGeno layers a "De-coding" layer on top of the Evo 2 foundation model:

1.  **Evo 2 (The Brain)**: A 7B parameter model that understands the grammar of all sequenced DNA on Earth.
2.  **Sparse Autoencoders (The Interpreter)**: 32,768 discrete concept nodes that translate internal AI math into human-interpretable biological signals.
3.  **Landmark Scouter**: Scans raw FASTQ data to find the precise coordinates where these concept nodes fire with the highest intensity.

---

## 🔬 Significance-First Workflow

PlatyGeno is organized as a professional, Ph.D.-grade discovery pipeline:

1. **Significance Scan**: `python validation/step1_discovery.py` — Identify all high-confidence landmarks.
2. **Novelty Validation**: `python validation/step2_local_blast.py` — identifies `novel_sequences.csv`.
3. **Structural Analysis**: `python validation/step3_fasta_prep.py` and `step4_alphafold_run.py`.

## ⚙️ Hardware Optimization

PlatyGeno is optimized for high-performance discovery. To resolve the "12-hour bottleneck" on large datasets, utilize the **Batched Inference** engine.

### Batch Size Guide (`--batch-size`)
Parallelizing your scan is the fastest way to get results. Match this setting to your GPU VRAM:

| Hardware | VRAM | Recommended Batch Size |
| :--- | :--- | :--- |
| **A100 / H100** | 80GB | `32` – `64` |
| **RTX 3090 / 4090** | 24GB | `8` – `16` |
| **RTX 3060 / 4070** | 12GB | `1` – `2` |

> [!TIP]
> **Out of Memory?** If you encounter an OOM error, simply lower the `--batch-size`.

---

## 🚀 Quick Start (Landmark Scan)

```bash
# 1. Install PlatyGeno
pip install platygeno

# 2. Run a Panoramic Significance Map
# This disables rarity filters to show you EVERY important landmark in the first 5000 reads.
platygeno --input sample.fastq --start 0 --limit 5000 --panoramic
```

---

## 🧩 The Scientific Dial: Tuning Significance

PlatyGeno uses the AI's "Excitement" as the primary scientific dial:

### 1. Signal Strength (`min_activation`)
*   **3.0 – 5.0**: "Significance Scouting." Ideal for mapping the general landscape of a sample.
*   **8.0 – 12.0**: "Landmark Identification." Targets high-confidence biological machinery.

### 2. Novelty Filter (`rel_freq_max`) - Optional
*   **Off (`1.0`)**: Standard mode. Shows all important genes (Known and Unknown).
*   **On (`0.001`)**: Novelty mode. Automatically subtracts common housekeeping genes to find "Dark Matter."

---

## 🧪 Use Case: Hunting for the Unknown

While PlatyGeno identifies all important genes, it is uniquely tuned for **Gromatic Dark Matter**:
*   **Reference-Free**: Identify significance in exotic metagenomes where no reference genomes exist.
*   **Structural Discovery**: Feed AI-flagged sequences directly into **AlphaFold** to discover never-before-seen 3D protein folds.

---

## 📚 API Reference

### `platygeno.discover_genes()`
| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `input_path` | `str` | *Req* | Path to sequence file. |
| `min_activation` | `float` | `5.0` | Minimum signal strength. |
| `rel_freq_max` | `float` | `1.0` | Rarity cap (1.0 = All significance). |
| `top_n` | `int` | `10` | Max number of landmarks to return. |

---

## 📜 Primary References

**1. PlatyGeno (This Package):**
```bibtex
@software{PlatyGeno2026,
  author = {Khoa Tu Tran},
  title = {PlatyGeno: Unsupervised Significance Mapping via Evo 2},
  url = {https://github.com/khoatran1995/PlatyGeno},
  year = {2026}
}
```

**2. Evo 2 Model:** Arc Institute. (2026). *Genome modeling and design across all domains of life with Evo 2*. *Nature*.
