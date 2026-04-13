import torch
import torch.nn as nn
from evo2 import Evo2
from huggingface_hub import hf_hub_download
import numpy as np

class SparseAutoencoder(nn.Module):
    def __init__(self, d_model=4096, d_hidden=32768):
        super().__init__()
        self.W_enc = nn.Parameter(torch.zeros(d_model, d_hidden))
        self.b_enc = nn.Parameter(torch.zeros(d_hidden))
        self.b_dec = nn.Parameter(torch.zeros(d_model))
        
    def encode(self, x, k=64):
        pre_act = (x - self.b_dec) @ self.W_enc + self.b_enc
        topk = torch.topk(torch.relu(pre_act), k, dim=-1)
        sparse_features = torch.zeros_like(pre_act)
        sparse_features.scatter_(-1, topk.indices, topk.values)
        return sparse_features

class PlatyGenoEngine:
    def __init__(self, model_name='evo2_7b'):
        self.device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
        print(f"🚀 Loading {model_name} into VRAM...")
        
        self.evo = Evo2(model_name)
        if self.device.startswith('cuda'):
            self.evo.model.half()
        self.evo.model.eval() 
        
        self.extracted_data = None
        self.evo.model.blocks[25].register_forward_hook(self._hook_fn)
        
        self.sae = self._setup_sae()
        print("✅ PlatyGeno Engine Ready.")

    def _hook_fn(self, module, input, output):
        data = output[0] if isinstance(output, tuple) else output
        self.extracted_data = data.detach()

    def _setup_sae(self):
        print("📡 Downloading/Loading SAE weights...")
        path = hf_hub_download(
            repo_id="Goodfire/Evo-2-Layer-26-Mixed",
            filename="sae-layer26-mixed-expansion_8-k_64.pt"
        )
        sae = SparseAutoencoder().to(self.device)
        state_dict = torch.load(path, map_location=self.device, weights_only=False)
        
        mapped_state = {
            'W_enc': state_dict['_orig_mod.W'],
            'b_enc': state_dict['_orig_mod.b_enc'],
            'b_dec': state_dict['_orig_mod.b_dec']
        }
        sae.load_state_dict(mapped_state)
        return sae.eval()

    def get_features(self, dna_strings):
        """
        High-Fidelity Feature Extraction.
        Processes sequences individually to avoid padding artifacts and maximize sensitivity.
        """
        if isinstance(dna_strings, str):
            dna_strings = [dna_strings]
            
        all_features = []
        for dna in dna_strings:
            try:
                tokens = self.evo.tokenizer.tokenize(dna)
                input_ids = torch.tensor([tokens], dtype=torch.long).to(self.device)
                
                with torch.no_grad():
                    self.extracted_data = None 
                    _ = self.evo.model(input_ids)
                    
                    if self.extracted_data is not None:
                        # Process token-level activations
                        token_features = self.sae.encode(self.extracted_data, k=64)
                        # Use Max-Pooling across the sequence to capture peak significance
                        read_features = torch.max(token_features, dim=1).values
                        all_features.append(read_features)
                    else:
                        all_features.append(torch.zeros((1, 32768), device=self.device))
            except torch.cuda.OutOfMemoryError:
                torch.cuda.empty_cache()
                all_features.append(torch.zeros((1, 32768), device=self.device))
                
        if not all_features:
            return None
            
        return torch.cat(all_features, dim=0)

    def get_token_features_deep(self, dna_string):
        tokens = self.evo.tokenizer.tokenize(dna_string)
        input_ids = torch.tensor([tokens], dtype=torch.long).to(self.device)
        with torch.no_grad():
            _ = self.evo.model(input_ids)
            if self.extracted_data is not None:
                token_embeddings = self.extracted_data.squeeze(0)
                return self.sae.encode(token_embeddings, k=64)
        return None
