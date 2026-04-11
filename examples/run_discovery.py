import torch
import _codecs
# Your Skeleton Key from the notebook
torch.serialization.add_safe_globals([_codecs.encode])

from platygeno.core import PlatyGenoEngine
from platygeno.evo_reader import read_evo_features

import os

# 1. Setup paths relative to the project root
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
fastq_file = os.path.join(base_dir, "data", "sample.fastq")
output_file = os.path.join(base_dir, "data", "novel_candidates.csv")

# 2. Start the Engine
engine = PlatyGenoEngine(model_name='evo2_7b')

# 3. Run Discovery on your FASTQ
if not os.path.exists(fastq_file):
    print(f"⚠️ Error: FASTQ file not found at {fastq_file}")
else:
    print(f"🔍 Scanning {fastq_file} for novel gene candidates...")
    report = read_evo_features(fastq_file, engine, limit=4000)

    # 4. Save results
    report.to_csv(output_file, index=False)
    print(f"🎉 Discovery Complete! Report saved to {output_file}")
