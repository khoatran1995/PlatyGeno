from platygeno.core import PlatyGenoEngine
from platygeno.evo_reader import read_evo_features
from platygeno.mapper import get_high_score_reads

# 1. Initialize with specific Evo version
engine = PlatyGenoEngine(sae_weights_path="sae.pt", model_version='evo-2-7b-base')

# 2. Run targeted range discovery
raw_report = read_evo_features(
    fasta_path="metagenome.fasta", 
    engine=engine,
    start_idx=0, 
    end_idx=4000, 
    score_thres=12.0
)

# 3. Extract high-importance signals
discovery = get_high_score_reads(raw_report, score_thres=20.0)
print(discovery.head())
