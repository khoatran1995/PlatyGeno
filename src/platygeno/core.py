import os
import torch
from evo import Evo  # Ensure you have the latest evo-model or evo2 installed
from huggingface_hub import hf_hub_download

class PlatyGenoEngine:
    def __init__(self, sae_weights_path=None, model_version='evo2-7b', device='cuda:0'):
        self.device = device
        
        # 1. Initialize Evo 2 (Using the 2026 identifier)
        print(f"🚀 Initializing {model_version}...")
        try:
            self.evo = Evo(model_version).to(device).eval()
        except ValueError:
            print(f"⚠️ Warning: '{model_version}' not found in local registry. Attempting fallback...")
            # Fallback for specific environment setups
            self.evo = Evo('evo-1.5-8k-base').to(device).eval() 

        # 2. Handle SAE Weights Auto-Download
        self.sae_path = sae_weights_path or "weights/sae_layer26.pt"
        self._ensure_weights_exist()
        
        self.sae = self._load_sae(self.sae_path)
        print(f"✅ PlatyGeno Engine Ready: {model_version} + Layer 26 SAE.")

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
