import os
import torch
from evo import Evo  # Ensure you have the latest evo-model or evo2 installed
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
                filename="sae_layer26.pt", # Verify exact filename on HF
                local_dir="weights"
            )
            self.sae_path = downloaded_path
            print(f"✅ Download complete: {self.sae_path}")

    def _load_sae(self, path):
        # Your specific SAE loading logic goes here
        # return torch.load(path, map_location=self.device)
        pass

    def get_features(self, dna_string):
        # Extraction logic...
        pass
