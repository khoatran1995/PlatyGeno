import torch
import _codecs
# Your Skeleton Key from the notebook
torch.serialization.add_safe_globals([_codecs.encode])

from platygeno.core import PlatyGenoEngine
from platygeno.evo_reader import read_evo_features

# 1. Start the Engine
engine = PlatyGenoEngine(model_name='evo2_7b')

# 2. Run Discovery on your FASTQ
# Ensure this path matches your RunPod workspace
fastq_file = "/workspace/PlatyGeno/data/sample.fastq"

print(f"🔍 Scanning {fastq_file} for novel gene candidates...")
report = read_evo_features(fastq_file, engine, limit=4000)

# 3. Save results
report.to_csv("novel_candidates.csv", index=False)
print("🎉 Discovery Complete! Report saved to novel_candidates.csv")
