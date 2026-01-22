import os
import pickle
import sys

import torch
from safetensors.torch import load_file, save_file

# -------------------------------
# Metadata encoding/decoding
# -------------------------------


def encode_metadata(obj):
    """
    Recursively encode Python objects into tensors:
    - torch.Tensor → leave as-is
    - dict → recursively encode
    - list/tuple → convert to dict {0: v0, 1: v1, ...} and encode recursively
    - other → pickle into uint8 tensor
    """
    if isinstance(obj, torch.Tensor):
        return obj
    elif isinstance(obj, dict):
        return {k: encode_metadata(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return {str(i): encode_metadata(v) for i, v in enumerate(obj)}
    else:
        data = pickle.dumps(obj)
        return torch.tensor(list(data), dtype=torch.uint8)


def decode_metadata(obj):
    """
    Recursively decode tensors back into Python objects.
    """
    if isinstance(obj, torch.Tensor):
        if obj.dtype == torch.uint8:
            data = bytes(obj.tolist())
            return pickle.loads(data)
        return obj
    elif isinstance(obj, dict):
        # Convert dicts with all digit keys back to lists
        if all(k.isdigit() for k in obj.keys()):
            return [decode_metadata(obj[k]) for k in sorted(obj.keys(), key=int)]
        else:
            return {k: decode_metadata(v) for k, v in obj.items()}
    else:
        return obj

# -------------------------------
# Flatten/unflatten for SafeTensors
# -------------------------------


def flatten_dict(d, parent_key='', sep='/'):
    items = {}
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.update(flatten_dict(v, new_key, sep=sep))
        else:
            items[new_key] = v
    return items


def unflatten_dict(d, sep='/'):
    result = {}
    for k, v in d.items():
        keys = k.split(sep)
        target = result
        for key in keys[:-1]:
            target = target.setdefault(key, {})
        target[keys[-1]] = v
    return result


# -------------------------------
# Save .pt as SafeTensors
# -------------------------------

if __name__ == "__main__":
    FILE_PATH = sys.argv[1]
    if FILE_PATH.endswith('.pt'):
        checkpoint = torch.load("model.pt", map_location="cpu")
        encoded = encode_metadata(checkpoint)
        flat = flatten_dict(encoded)
        save_file(flat, os.path.join(os.path.dirname(sys.argv[1]), "model.safetensors"))
        print("Saved restorable SafeTensors file!")
    else:
        loaded_flat = load_file("model_restorable.safetensors")
        loaded_nested = unflatten_dict(loaded_flat)
        restored_checkpoint = decode_metadata(loaded_nested)
        torch.save(restored_checkpoint, os.path.join(os.path.dirname(sys.argv[1]), "model.pt"))
        print("Saved restored checkpoint as model_restored.pt!")
