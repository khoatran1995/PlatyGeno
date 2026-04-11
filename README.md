# PlatyGeno 🧬
PlatyGeno is an unsupervised gene discovery package designed to interpret the **Evo 2 genomic foundation model** using Sparse Autoencoders (SAEs).

It identifies functional biological motifs (promoters, coding regions, etc.) by mapping high-activation internal features back to the original DNA sequences.

---

## 🚀 Installation — RunPod (Recommended)

### Step 1 — Clone the Repo
```bash
git clone https://github.com/khoatran1995/PlatyGeno.git
cd /workspace/PlatyGeno
```

### Step 2 — Environment Setup
```bash
pip install -e .
```

### Step 3 — Run the Discovery
You can choose between an automated pipeline or a step-by-step custom workflow:

#### 1. All-in-One Pipeline
The automated "1-for-all" discovery script. Uses the high-level API.
```bash
python examples/all_in_one_discovery.py --start 0 --end 4000
```

#### 2. Step-by-Step Workflow
A manual, modular look at every phase of the discovery logic.
```bash
python examples/step_by_step_discovery.py --start 0 --end 4000
```

---

## 💻 Library Usage (One-Line API)

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
| `method` | **Best Snippet** or **Assembled Contig**. |
| `feature_id` | The biological motif identified by the model. |
| `read_id` | Source read or "Consensus" (if assembled). |
| `activation` | Strength of the biological signal. |
| `length` | Result length in base-pairs. |
| `sequence` | The final DNA sequence result. |

---

## 🛠 Dependencies
This project requires access to:
- **Evo 2 weights** (Arc Institute)
- **SAE weights** (Goodfire)

## Acknowledgements
Developed by Khoa Tu Tran. Built upon the Evo 2 foundational work by the Arc Institute and SAE architectures by the Goodfire team.
