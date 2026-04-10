from platygeno.core import PlatyGenoEngine
from platygeno.evo_reader import read_evo_features
from platygeno.mapper import get_high_score_reads

# 1. Initialize Engine
# SET sae_weights_path=None to trigger the automatic download logic we wrote!
# Use 'evo2-7b' (standard ID) or ensure your core.py handles the fallback.
engine = PlatyGenoEngine(
    sae_weights_path=None, 
    model_version='evo2-7b' 
)

# 2. Run targeted range discovery
# Update fasta_path to point to your 'data/' folder
raw_report = read_evo_features(
    fasta_path="data/sample.fastq",  # Change .fasta to .fastq
    engine=engine,
    start_idx=0, 
    end_idx=4000, 
    score_thres=12.0
)

# 3. Extract high-importance signals
discovery = get_high_score_reads(raw_report, score_thres=20.0)

# 4. Final output
if not discovery.empty:
    print(discovery.head())
    discovery.to_csv("discovery_results.csv", index=False)
    print("✅ Results saved to discovery_results.csv")
else:
    print("❌ No features found above the threshold.")
