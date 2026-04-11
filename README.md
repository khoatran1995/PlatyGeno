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

This compiles the CUDA kernel for FlashAttention. It takes ~5–10 minutes:

```bash
pip install flash-attn --no-build-isolation
```

> ℹ️ The `WARNING: Running pip as 'root'` message is expected on RunPod — it is safe to ignore.

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

```bash
# Example: download from a public URL
wget -O data/sample.fastq "YOUR_FILE_URL"
```

### Step 6 — Run the Discovery

```bash
python examples/gene_discovery_sample.py
```

> ⚠️ Do **not** run as `examples/gene_discovery_sample.py` — you will get a "Permission denied" error.
> Always prefix with `python`.

---

## 📁 Output Files

Results are saved to the `data/` directory:

```
data/
├── phase1_feature_report.csv    ← Full feature activation report (all reads)
└── gene_snippets_top10.csv      ← Top 10 discovered gene snippets
```

The final `gene_snippets_top10.csv` contains:

| Column | Description |
|--------|-------------|
| `rank` | Discovery rank (1 = strongest signal) |
| `feature_id` | SAE feature index |
| `read_id` | Source read from the input file |
| `phase1_activation` | Mean-pooled activation from Phase 1 |
| `peak_activation` | Precise token-level peak activation |
| `peak_bp_index` | Base-pair position of the peak in the read |
| `gene_snippet_60bp` | The extracted 60bp DNA sequence |

---

## 🐛 Troubleshooting

| Error | Cause | Fix |
|-------|-------|-----|
| `No module named 'flash_attn_2_cuda'` | FlashAttention not installed | Run `pip install flash-attn --no-build-isolation` |
| `ModuleNotFoundError: evo2` | Evo 2 not installed | Run `pip install -e .` |
| `FileNotFoundError: sample.fastq` | Input file missing | Upload your FASTQ to `data/sample.fastq` |
| `CUDA out of memory` | GPU too small | Use A100 80GB or reduce `SCAN_LIMIT` in the script |
| `Permission denied` when running script | Missing `python` prefix | Use `python examples/gene_discovery_sample.py` |
| HuggingFace rate limit warning | No HF token set | Run `export HF_TOKEN=your_token` before the script |

---

## Acknowledgements

This tool is built upon the **Evo 2 genomic foundation model** developed by the Arc Institute.
We are grateful to the Arc Institute and the Goodfire team for providing the model weights and
sparse autoencoder (SAE) architectures that make PlatyGeno possible.

This project was developed with the assistance of AI tools for code optimization and debugging,
with all biological logic and architectural design directed by the author.
