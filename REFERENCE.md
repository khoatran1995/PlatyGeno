# PlatyGeno | API Reference (v1.0.0)

This document provides the technical reference for the `platygeno` Python package. Use these functions to build custom discovery pipelines or integrate PlatyGeno into Jupyter Notebooks.

---

## 🚀 Master Workflow

### `platygeno.discover_genes()`
The primary end-to-end mission controller for biological significance mapping.

**Parameters:**
- `input_path` (str): Path to the FASTQ/FASTA file.
- `engine` (PlatyGenoEngine, optional): An existing engine instance to prevent re-loading weights.
- `scan_start` (int): Starting read index. Default: `0`.
- `scan_end` (int): Ending read index. Set to `None` for full file scan.
- `top_n` (int): Number of features to report. Set to `-1` for all significant hits.
- `min_activation` (float): Sensitivity threshold. Default: `1.0`.
- `rel_freq_max` (float): Rarity cap. Set to `1.0` for panoramic discovery.
- `batch_size` (int): GPU batch size (Higher = Faster).

**Returns:**
- `pd.DataFrame`: A comprehensive report containing significance peaks, species origins, and sequence snippets.

---

## 🧬 Core Engine

### `platygeno.PlatyGenoEngine()`
Initializes the Evo 2 foundation model and the Sparse Autoencoder (SAE) weights.

**Key Methods:**
- `.get_features(dna_strings)`: Extracts biological signals for a batch of sequences.
- `.get_token_features_deep(dna_string)`: Performs token-level scanning for precise boundary detection.

---

## 📡 Data Extraction

### `platygeno.read_evo_features()`
The high-speed reader that manages GPU memory and forward passes.

**Parameters:**
- `file_path` (str): Path to DNA sequence data.
- `engine` (PlatyGenoEngine): Initialized engine instance.
- `start / stop` (int): Coordinate range in the file.
- `batch_size` (int): Number of reads processed per GPU pass.

---

## 📊 Significance Mapping

### `platygeno.find_significant_landmarks()`
Identifies biological points of interest based on activation strength and population statistics.

**Parameters:**
- `df` (pd.DataFrame): The raw report from `read_evo_features`.
- `min_activation` (float): The "volume" threshold for a biological signal.
- `total_population` (int): Total reads scanned (used for rarity math).

### `platygeno.annotate_with_biology()`
Cross-references discovered Feature IDs with the internal `layer26_features.csv` dictionary.

---

## 📈 Reporting

### `platygeno.generate_html_report()`
Generates the premium, dark-mode Discovery Dashboard.

**Parameters:**
- `results_df` (pd.DataFrame): The final annotated discovery dataframe.
- `output_path` (str): Destination for the HTML file.
