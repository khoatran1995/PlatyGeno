# Copyright 2026 Khoa Tu Tran
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pandas as pd
from Bio import SeqIO
from platygeno import PlatyGenoEngine

# --- Configuration ---
FASTQ_FILE = "data/sample.fastq"
SAE_PATH = "weights/sae_layer26.pt"
OUTPUT_FILE = "discovery_report.csv"
# ---------------------

engine = PlatyGenoEngine(SAE_PATH)
findings = []

print(f"Analyzing {FASTQ_FILE}...")

for record in SeqIO.parse(FASTQ_FILE, "fastq"):
    activations = engine.get_features(str(record.seq))
    
    # Identify top active features
    indices = activations.nonzero(as_tuple=True)[0]
    for idx in indices:
        val = activations[idx].item()
        if val > 10.0:  # Significance threshold
            findings.append({
                "read_id": record.id,
                "feature_id": idx.item(),
                "activation": val
            })

df = pd.DataFrame(findings)
df.to_csv(OUTPUT_FILE, index=False)
print(f"Success. Results saved to {OUTPUT_FILE}")
