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
        # Your winning math from the notebook
        pre_act = (x - self.b_dec) @ self.W_enc + self.b_enc
        topk = torch.topk(torch.relu(pre_act), k, dim=-1)
        sparse_features = torch.zeros_like(pre_act)
        sparse_features.scatter_(-1, topk.indices, topk.values)
        return sparse_features

class PlatyGenoEngine:
    def __init__(self, model_name='evo2_7b'):
        self.device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
        print(f"🚀 Loading {model_name} into VRAM...")
        
        # Load Evo 2 (Using the wrapper class)
        # Optimized for 24GB VRAM (3090/4090)
        self.evo = Evo2(model_name)
        if self.device.startswith('cuda'):
            self.evo.model.half() # Move to half-precision to save 50% VRAM
        self.evo.model.eval() 
        
        # Setup Hook for Layer 26 (Index 25)
        self.extracted_data = None
        self.evo.model.blocks[25].register_forward_hook(self._hook_fn)
        
        # Load SAE
        self.sae = self._setup_sae()
        print("✅ PlatyGeno Engine Ready.")

    def _hook_fn(self, module, input, output):
        # Logic from your notebook to capture hidden states
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
        
        # The "Ghost Key" mapping from your notebook
        mapped_state = {
            'W_enc': state_dict['_orig_mod.W'],
            'b_enc': state_dict['_orig_mod.b_enc'],
            'b_dec': state_dict['_orig_mod.b_dec']
        }
        sae.load_state_dict(mapped_state)
        return sae.eval()

    def get_features(self, dna_strings):
        """
        Batched Feature Extraction with Recursive OOM Protection.
        If a batch fails due to VRAM, it splits and retries.
        """
        if not isinstance(dna_strings, (list, tuple)):
            dna_strings = [dna_strings]

        if not dna_strings:
            return None

        # 1. Batched Tokenization
        try:
            token_list = [torch.tensor(self.evo.tokenizer.tokenize(s), dtype=torch.long) for s in dna_strings]
            input_ids = torch.nn.utils.rnn.pad_sequence(token_list, batch_first=True, padding_value=0).to(self.device)
            
            # 2. Forward pass to trigger hook
            with torch.no_grad():
                self.extracted_data = None 
                _ = self.evo.model(input_ids)
                
                if self.extracted_data is not None:
                    mean_emb = torch.mean(self.extracted_data, dim=1)
                    features = self.sae.encode(mean_emb, k=64)
                    return features
                    
        except torch.cuda.OutOfMemoryError:
            torch.cuda.empty_cache()
            batch_size = len(dna_strings)
            
            # If we OOM with a single sequence, skip it and notify the user
            if batch_size <= 1:
                seq_len = len(dna_strings[0]) if dna_strings else "Unknown"
                print(f"⚠️ SKIPPING: Sequence too long (Length: {seq_len}) for current GPU VRAM.")
                # Return a zero-tensor of the correct hidden dimension (32768) to maintain alignment
                return torch.zeros((1, 32768), device=self.device)
            
            # Otherwise, split the batch and try again
            mid = batch_size // 2
            print(f"⚠️ GPU OOM: Splitting batch ({batch_size} -> {mid}) and retrying...")
            
            part1 = self.get_features(dna_strings[:mid])
            part2 = self.get_features(dna_strings[mid:])
            
            # Since skipped single sequences now return zeros instead of None, 
            # we can safely concatenate the results back together.
            return torch.cat([part1, part2], dim=0)

        return None

    def get_token_features_deep(self, dna_string):
        """
        Phase 3: Deep Token-Aware Scan.
        Returns SAE features for EVERY base pair (token) in the sequence.
        Used only on 'Winner' reads to find the precise gene boundaries.
        Returns shape: [Seq_Len, d_hidden]
        """
        tokens = self.evo.tokenizer.tokenize(dna_string)
        input_ids = torch.tensor([tokens], dtype=torch.long).to(self.device)

        with torch.no_grad():
            _ = self.evo.model(input_ids)
            if self.extracted_data is not None:
                # Squeeze to get [Seq_Len, Hidden_Dim]
                token_embeddings = self.extracted_data.squeeze(0)
                # Apply SAE to every single base token
                return self.sae.encode(token_embeddings, k=64)
        return None
