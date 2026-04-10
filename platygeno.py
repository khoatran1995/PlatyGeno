# Copyright 2026 Khoa Tu Tran
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

import torch
import torch.nn as nn
from evo import Evo

class SparseAutoencoder(nn.Module):
    """JumpReLU Sparse Autoencoder for Layer 26 activations."""
    def __init__(self, d_model=4096, d_sae=32768):
        super().__init__()
        self.W_enc = nn.Parameter(torch.zeros(d_model, d_sae))
        self.b_enc = nn.Parameter(torch.zeros(d_sae))
        self.b_dec = nn.Parameter(torch.zeros(d_model))

    def encode(self, x, threshold=0.0):
        # x: (batch, d_model)
        # Centering and encoding
        pre_act = (x - self.b_dec) @ self.W_enc + self.b_enc
        # JumpReLU activation
        return torch.where(pre_act > threshold, pre_act, torch.zeros_like(pre_act))

class PlatyGenoEngine:
    def __init__(self, sae_weights_path, device='cuda:0'):
        self.device = device
        self.layer = 26
        
        # Initialize Evo 2 (7B)
        print("Loading Evo 2...")
        self.evo = Evo('evo-2-7b-base').to(device).eval()
        
        # Initialize and Load SAE
        self.sae = SparseAutoencoder().to(device).eval()
        state_dict = torch.load(sae_weights_path, map_location=device)
        
        # Weight mapping (Assumes Goodfire/Arc naming conventions)
        self.sae.load_state_dict({
            'W_enc': state_dict['_orig_mod.W'],
            'b_enc': state_dict['_orig_mod.b_enc'],
            'b_dec': state_dict['_orig_mod.b_dec']
        })
        print("PlatyGeno Engine Ready.")

    def get_features(self, dna_sequence):
        with torch.no_grad():
            input_ids = self.evo.tokenizer.encode(dna_sequence, return_tensors='pt').to(self.device)
            _, hidden_states = self.evo(input_ids, output_hidden_states=True)
            
            # Mean pool the activations of Layer 26
            latent = hidden_states[self.layer].mean(dim=1) 
            features = self.sae.encode(latent)
            return features.squeeze()
