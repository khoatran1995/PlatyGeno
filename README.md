<table width="100%">
  <tr>
    <td width="200"><img src="https://raw.githubusercontent.com/khoatran1995/PlatyGeno/main/icon/PlatyGeno_v3.jpg" width="200" alt="PlatyGeno Icon"></td>
    <td align="left"><h1>Unsupervised Biological Significance Mapping via <br> Evo 2 & Sparse Autoencoders</h1></td>
  </tr>
</table>

[![PyPI version](https://img.shields.io/pypi/v/platygeno.svg)](https://pypi.org/project/platygeno/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19581440.svg)](https://doi.org/10.5281/zenodo.19581440)

PlatyGeno identifies **genomic landmarks** directly from raw sequence data. By leveraging the **Evo 2 foundation model**, it isolates biologically significant DNA structures (promoters, coding sequences, precise motifs) based purely on AI confidence—**without requiring labels, databases, or BLAST.**

## 🧪 Technical Validation (IBD-MDB)
PlatyGeno v1.0.3 has been validated using the clinical **[IBD Metagenomic Database](https://ibdmdb.org/)** dataset. For complete statistical data and methodology, see the [**PlatyGeno Technical Audit**](./PlatyGeno%20Technical%20Audit.md) and the [**PlatyGeno Technical Supplemental**](./PlatyGeno%20Technical%20Supplemental.md).

### Summary of Results:
*   **Novel Genomic Landmarks**: Identified **Feature 7393**, a 101bp sequence with no prior database matches and a high-confidence structural model (AlphaFold2 best prediction; $pLDDT \sim 80$).

    <p align="center">
      <img src="validation/sample_validation_analysis/feature7393_Alphafold2_best_structure.png" width="450" alt="Feature 7393 Structure">
      <br><i>Representative best structural prediction (AlphaFold2) of the novel Feature 7393 discovered autonomously by PlatyGeno.</i>
    </p>
*   **Statistical Correlation**: Verified a Pearson correlation of **$r=0.84$** ($p < 10^{-50}$) between sequence length and match significance.
*   **Resolution Gain**: Consensus assembly provided a **$10^{38}$ increase** in E-value confidence over isolated 60bp fragments.
*   **Taxonomic Profile**: 72% of high-activation discoveries successfully cross-validated with target gut microbiota.

---

## 🏗️ Technical Foundation

PlatyGeno operates as a **Reference-Free Microscope**, detecting the "Signal" of life directly from genomic grammar.

### 🔭 The Discovery Core
- **AI-Native Interpretation**: We use a **Sparse Autoencoder (SAE)** to translate the complex DNA "grammar" understood by **Evo 2** into 32,768 human-interpretable biological concepts (e.g., *promoters*, *viral motifs*).
- **Peak Pinpointing (Layer 26)**: The engine intercepts signals at Layer 26 to identify the exact coordinate where a biological feature fires with the highest intensity.
- **Dual-Mode Discovery**: Preserves both narrow **Precision Snippets** (separatedly high-interest DNA clips) and **Consensus Assemblies** (overlapping sequences from multiple reads of the same feature pieced together ).

> [!IMPORTANT]
> **Performance Highlight**: While both modes are preserved in discovery, validation benchmarks confirm that **Consensus Assembly** yields statistically superior significance (E-values) and cleaner taxonomic resolution.

👉 **For a full hierarchical deep-dive into the methodology and validation trail, see [Technical Architecture](docs/ARCHITECTURE.md).**

---

## ⚙️ Setup & Installation

### ⚙️ Installation & Quick Start

PlatyGeno requires a CUDA-enabled GPU (RTX 3090, 4090, A100, or H100).

```bash
# 1. Install the core package
pip install platygeno

# 2. Install high-performance GPU kernels (Mandatory for speed)
pip install ninja # for faster installation of flash-attn
pip install flash-attn --no-build-isolation

# 3. Verify & Run Discovery (on the validation sample)
platygeno --input data/sample.fastq --limit 5000 --threshold 5.0
```

### 🚀 Quick Start for GitHub Clones
```bash
# 1. Clone & Enter
git clone https://github.com/khoatran1995/PlatyGeno.git
cd PlatyGeno

# 2. Install High-Performance Kernels & editable package
pip install flash-attn --no-build-isolation
pip install -e .

# 3. Trigger Discovery
platygeno --input data/sample.fastq --limit 5000
```

### 📚 Documentation
**[API Reference](docs/DOCUMENTATION.md)**: Details on Evo 2 integrations and technical Python parameters.

---

## 🚀 Usage & API Reference

### 🚀 Advanced Python Discovery
Researchers can integrate the engine into custom discovery pipelines:

```python
import platygeno

# Advanced Discovery: Tuning parameters for clinical audits
results = platygeno.discover_genes(
    input_path="data/sample.fastq",
    scan_end=5000,
    min_activation=8.0,      # High-confidence threshold
    batch_size=32            # GPU-optimized batching
)

# View discovered biological features
print(results[['feature_id', 'feature_name', 'activation', 'sequence']])
```

### `platygeno.discover_genes()` Reference
| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `input_path` | `str` | *Req* | Path to sequence file. |
| `min_activation` | `float` | `5.0` | Minimum signal strength. |
| `rel_freq_max` | `float` | `1.0` | Rarity cap (1.0 = All significance). |
| `scan_end` | `int` | `None` | Last read index (**None for end of file**). |
| `top_n` | `int` | `-1` | Max features to return (**-1 for ALL**). |

---

## ⚡ Performance
| Mode | Engine Implementation | Runtime (20k Reads) | Discovery Speed |
| :--- | :--- | :--- | :--- |
| **v1.0.3** | **Batched Mean-Pooling** | **~4.8 Minutes** | **🚀 100% (High Speed)** |

---

## 📜 References

## ⚠️ Technical Limitations
*   **Pre-training Bias**: Sensitivity depends on the **Evo 2** pre-training corpus.
*   **SAE Bottleneck**: Discrete compression may miss extremely subtle biological nuances.
*   **Validation Requirement**: High significance is a "Beacon," not final functional proof.

---

## 📜 References
```bibtex
@software{PlatyGeno2026,
  author = {Khoa Tu Tran},
  title = {PlatyGeno: Unsupervised Significance Mapping via Evo 2},
  url = {https://github.com/khoatran1995/PlatyGeno},
  doi = {10.5281/zenodo.19581440},
  year = {2026}
}
```
*Thanks to Together AI (Evo 2) and Goodfire AI (SAE interpretability) and the IBD Metagenomic Database (IBD-MDB). Please cite the relevant references when using PlatyGeno.*
