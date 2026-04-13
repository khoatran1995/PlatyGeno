<table width="100%">
  <tr>
    <td width="200"><img src="icon/PlatyGeno%20-%20small.png" width="200" alt="PlatyGeno Icon"></td>
    <td align="left"><h1>Unsupervised Biological Significance Mapping via <br> Evo 2 & Sparse Autoencoders</h1></td>
  </tr>
</table>

[![PyPI version](https://img.shields.io/pypi/v/platygeno.svg)](https://pypi.org/project/platygeno/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

PlatyGeno is a professional Python package for identifying **genomic landmarks** directly from raw sequence data. By leveraging the **Evo 2 foundation model**, it identifies biologically significant DNA structures (promoters, coding sequences, precise motifs) based purely on AI confidence—**without requiring labels, databases, or BLAST.**

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
pip install ninja # for faster installation of flash-attn
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

## 📚 Documentation & Reference
*   **[Technical API Reference](REFERENCE.md)**: Detailed documentation for every function in the `platygeno` core.
*   **[Architecture Guide](ARCHITECTURE.md)**: Deep dive into Evo 2, Sparse Autoencoders, and Max-Pooling theory.

---

## 🐍 Python API Quickstart
Researchers can use PlatyGeno directly in Python/Jupyter for custom discovery workflows:

```python
import platygeno

# End-to-End Discovery
results = platygeno.discover_genes(
    input_path="data/sample.fastq",
    scan_end=1000,          # Scan first 1000 reads
    min_activation=1.0      # High sensitivity
)

# Access the discovered sequences
print(results[['feature_id', 'feature_name', 'sequence']])
```

---

## 🚀 Step-by-Step Discovery (Ph.D. Suite)

PlatyGeno is now organized as a unified, Ph.D.-grade discovery workflow:

1. **One-Touch Discovery**: `python validation/discovery_pipeline.py --input sample.fastq` — Performs both significance scanning and automated BLAST validation (via `validation/step2_blast.py`).
2. **AI-Aware Validation**: The engine automatically labels known features (Coding Regions, Alpha Helices) and prioritizes unknown "Dark Matter" for validation.
3. **Recursive OOM Guard**: Automatically scales batch sizes to fit your GPU VRAM, ensuring large files don't crash the discovery process.

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

# 2. Run a Significance Map
# This will show you EVERY significant landmark in the first 5000 reads.
platygeno --input sample.fastq --limit 5000

# To Scan the entire file (Start to End):
platygeno --input sample.fastq --limit -1
```

---

## 🧩 The Scientific Dial: Tuning Significance

PlatyGeno uses the AI's "Excitement" as the primary scientific dial:

### 1. Signal Strength (`min_activation`)
*   **3.0 – 5.0**: "Significance Scouting." Ideal for mapping the general landscape of a sample.
*   **8.0 – 12.0**: "Landmark Identification." Targets high-confidence biological machinery.

### 2. Novelty Filter (`--rarity-only`) - Optional
*   **Default (Off)**: Standard mode (Panoramic). Shows all important genes (Known and Unknown).
*   **On (`--rarity-only`)**: Novelty mode. Automatically subtracts common housekeeping genes to find "Dark Matter."

### 3. Discovery Breadth (`--top-n`)
*   **Default (`-1`)**: Ph.D. Survey Mode. Returns **every significant landmark** found in the sample (Unlimited).
*   **Targeted (`10-25`)**: Precision Mode. Focuses only on the strongest outliers.

## 🧪 Core Methodology (The "Easy" Version)

PlatyGeno’s "Golden Configuration" is built on two simple rules:

### 1. Mean-Pooling (The "Big Picture")
Instead of getting distracted by every tiny detail in a DNA sequence, the AI looks at the **whole read at once**. It summarizes the sequence so it can find the main biological theme without getting lost in the noise.

### 2. Zero-Gate Discovery (The "No-Filter" Mode)
Most tools only show you the "Top 10" most interesting things they find. **Zero-Gate** means we don't filter anything out—if the AI sees it, you see it. This ensures that even the rarest, hidden genes find their way into your report.

---

### 🧬 Technical Performance Highlights (v1.0)
*   **Mechanism**: **Mean-Pooling** (Iterative sequence averaging).
*   **Diversity**: **Zero-Gate Discovery** (Captures ALL active biological signals).
*   **Performance**: Optimized for the **98-hit** Ph.D. benchmark discovery.

### 4. Strategic Subtraction (`--exclude`)
*   **Usage**: `--exclude 212,16509`
*   **Purpose**: "Mutes" features you have already identified as known biology. This forces the engine to look deeper and surface the next layer of genomic candidates.

## 🧪 Benchmark Dataset: Gut Metagenome (IBD-MDB)
PlatyGeno includes a high-density clinical benchmark for testing novelty discovery in complex human samples:

*   **Origin**: Chronic Inflammatory Bowel Disease (IBD) Metagenomic Database.
*   **Role**: Validating the engine's ability to identify autonomous biological landmarks in high-complexity clinical metagenomes.
*   **Local Data**: Benchmark reads are provided in the `data/` directory for Ph.D. reproducibility.

---

## 🧪 Use Case: Hunting for the Unknown

While PlatyGeno identifies all important genes, it is uniquely tuned for **Genomic Dark Matter**:
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
| `scan_end` | `int` | `None` | Last read index (**None for end of file**). |
| `top_n` | `int` | `-1` | Max features to return (**-1 for ALL, Default**). |

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
