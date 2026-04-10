import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

def analyze_reads(fasta_path, engine, limit=100, layer=26):
    """Function 1: Run reads through Evo 2 and return Feature Report."""
    results = []
    
    for i, record in enumerate(tqdm(SeqIO.parse(fasta_path, "fasta"))):
        if i >= limit: break
        
        # Get activations from core engine
        activations = engine.process(str(record.seq), layer=layer)
        
        # Find features above zero/threshold
        indices = activations.nonzero(as_tuple=True)[0]
        for idx in indices:
            results.append({
                "read_id": record.id,
                "seq": str(record.seq),
                "feature_id": idx.item(),
                "strength": activations[idx].item()
            })
            
    return pd.DataFrame(results)
