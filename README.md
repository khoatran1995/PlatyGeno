# PlatyGeno 🧬

PlatyGeno is an interpretability tool for the **Evo 2 genomic foundation model**.
It uses Sparse Autoencoders (SAEs) to perform unsupervised, token-level discovery of
functional gene candidates directly from raw FASTA/FASTQ sequencing reads.

---

## 🔬 How It Works

PlatyGeno runs a 3-phase discovery pipeline:

| Phase | Description |
|-------|-------------|
| **Phase 1 — Discovery Scan** | Streams reads through Evo 2 + SAE. Outputs a feature activation report (`read_id`, `feature_id`, `activation`). |
| **Phase 2 — Rare Signal Filtering** | Applies statistical filters to isolate rare, high-activation biological signals. |
| **Phase 3 — Gene Snippet Extraction** | Re-runs a deep token-level scan on the winning reads. Extracts the exact 60bp nucleotide window responsible for the highest activation. |

---

## 🛠 External Dependencies

These are **not included** in this repository and must be installed separately:

- **Evo 2 Foundation Model** (Arc Institute): https://huggingface.co/arcinstitute
- **Sparse Autoencoder Weights** (Goodfire): https://huggingface.co/Goodfire

Both are downloaded automatically on first run.

---

## 🚀 Installation — RunPod (Recommended)

> PlatyGeno requires a GPU with at least **40GB VRAM**.
> Recommended: **NVIDIA A100 80GB** on RunPod with the `RunPod PyTorch` template.

### Step 1 — Create a RunPod Instance

- **GPU**: A100 80GB (minimum)
- **Template**: RunPod PyTorch (Python 3.11 + CUDA pre-installed)
- **Disk**: 100GB container disk

### Step 2 — Clone this Repository

Open the RunPod terminal and run:

```bash
cd /workspace
git clone https://github.com/khoatran1995/PlatyGeno.git
cd PlatyGeno
```

### Step 3 — Install FlashAttention (Required for Evo 2)

This compiles the CUDA kernel for FlashAttention:

```bash
pip install flash-attn --no-build-isolation
```

### Step 4 — Install PlatyGeno and All Dependencies

```bash
pip install -e .
```

This installs: `torch`, `numpy`, `biopython`, `pandas`, `tqdm`, `huggingface_hub`, and `evo2`.

### Step 5 — Add Your Input File

Upload your `.fastq` or `.fasta` file to the pod. The default expected path is:

```
/workspace/PlatyGeno/data/sample.fastq
```

You can upload via the **JupyterLab file browser** (drag & drop), or use:

```
# Example: download from a public URL
wget -O data/sample.fastq "YOUR_FILE_URL"
```

---

## 🚀 Installation — RunPod (Recommended)

...

### Step 6 — Run the Discovery

You can choose between a detailed step-by-step example or a simplified automated pipeline:

#### 1. Detailed Manual Workflow (Learn Mode)
Shows exactly how Phase 1, 2, and 3 work under the hood.
```bash
python examples/gene_discovery_sample.py --start 0 --end 4000
```

#### 2. Automated Pipeline (Production Mode)
Uses the high-level `platygeno.discover_genes()` API for ease of use.
```bash
python examples/gene_discovery_pipeline.py --start 10000 --end 20000 --threshold 15.0
```

---

## 💻 Library Usage (One-Line API)

You can use PlatyGeno in your own Python scripts with a single function call:

```python
import platygeno

# Returns a DataFrame and optionally saves to CSV
results = platygeno.discover_genes(
    input_path="data/sample.fastq",
    scan_start=0,
    scan_end=5000,
    min_activation=10.0,
    output_path="results.csv"
)
```

---

## 📁 Output Files

Results are saved to `data/gene_discovery_results.csv`:

| Column | Description |
|--------|-------------|
| `method` | **Best Snippet** (single read) or **Assembled Contig** (merged consensus). |
| `feature_id` | The biological motif identified by the model. |
| `read_id` | The source read or "Consensus" (if assembled). |
| `activation` | The strength of the biological signal. |
| `length` | Result length in base-pairs. |
| `sequence` | The final DNA sequence result. |


## Acknowledgements

This tool is built upon the **Evo 2 genomic foundation model** developed by the Arc Institute.
We are grateful to the Arc Institute and the Goodfire team for providing the model weights and
sparse autoencoder (SAE) architectures that make PlatyGeno possible.

This project was developed with the assistance of AI tools for code optimization and debugging,
with all biological logic and architectural design directed by the author.
