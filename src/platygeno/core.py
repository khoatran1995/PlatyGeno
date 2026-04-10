import os
import torch
from evo2 import Evo2  # Ensure you have the latest evo-model or evo2 installed
from huggingface_hub import hf_hub_download

class PlatyGenoEngine:
    def __init__(self, sae_weights_path=None, model_name='evo2_7b'):
        # Evo2 handles device placement automatically. 
        # It will use CUDA if a GPU is detected.
        print(f"🚀 Initializing {model_name}...")
        
        # Load Evo 2 (The correct ID is 'evo2_7b')
        self.evo = Evo2(model_name) 
        
        # Handle SAE Weights
        self.sae_path = sae_weights_path or "weights/sae_layer26.pt"
        self._ensure_weights_exist()
        
        # self.sae = self._load_sae(self.sae_path)
        print(f"✅ PlatyGeno Engine Ready: {model_name} loaded.")


    def _ensure_weights_exist(self):
        """Checks if SAE weights exist locally; if not, downloads from Hugging Face."""
        if not os.path.exists(self.sae_path):
            print(f"📥 Weights not found at {self.sae_path}. Downloading from Goodfire/HuggingFace...")
            os.makedirs(os.path.dirname(self.sae_path), exist_ok=True)
            
            # This fetches the official Goodfire Layer 26 SAE
            downloaded_path = hf_hub_download(
                repo_id="Goodfire/Evo-2-Layer-26-Mixed",
                filename="sae-layer26-mixed-expansion_8-k_64.pt", # Verify exact filename on HF
                local_dir="weights"
            )
            self.sae_path = downloaded_path
            print(f"✅ Download complete: {self.sae_path}")


    def _load_sae(self, path):
        """Loads the SAE weights into memory."""
        print(f"📂 Loading SAE weights from {path}...")
        # Most Goodfire SAEs are saved as standard torch state dicts
        # We load them and move them to the GPU
        sae_data = torch.load(path, map_location=self.device)
        return sae_data

    def get_features(self, dna_string):
        """Processes a DNA string and returns SAE feature activations."""
        # 1. Convert DNA string to tokens
        # Evo 2 handles strings directly or via its internal tokenizer
        input_ids = self.evo.tokenize(dna_string).to(self.device)
        
        # 2. Run through Evo 2 and stop at Layer 26
        with torch.no_grad():
            # We get the 'hidden states' (activations) from layer 26
            outputs = self.evo(input_ids, return_embeddings=True, layer_idx=26)
            activations = outputs # This depends on the exact return shape of Evo2
            
            # 3. Pass activations through the SAE
            # If your SAE is a simple matrix multiplication (Encoder):
            # feature_acts = torch.relu(activations @ self.sae['encoder.weight'].T + self.sae['encoder.bias'])
            
            # For now, let's return a dummy tensor if you're still debugging 
            # to make sure the reader works:
            return torch.zeros(len(dna_string)) # Replace with real SAE logic
