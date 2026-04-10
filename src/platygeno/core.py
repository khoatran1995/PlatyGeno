import torch
from evo import Evo

class EvoModelWrapper:
    def __init__(self, model_name='evo-2-7b-base', device='cuda'):
        print(f"Initializing {model_name}...")
        self.device = device
        self.model = Evo(model_name).to(device).eval()
        self.tokenizer = self.model.tokenizer

    def get_layer_activations(self, sequence, layer=26):
        inputs = self.tokenizer.encode(sequence, return_tensors='pt').to(self.device)
        with torch.no_grad():
            _, hidden_states = self.model(inputs, output_hidden_states=True)
        return hidden_states[layer]
