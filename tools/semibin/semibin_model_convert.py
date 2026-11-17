"""Convert between the PyTorch checkpoints used by SemiBin and safetensors

SemiBin trains its models as PyTorch checkpoints (a ``dict`` with the entries
``model_name``, ``model_state_dict``, and ``params``), while Galaxy stores and
displays them as ``safetensors`` files with the same information in the
safetensors metadata section.

Usage::

    semibin_model_convert.py pt2safetensors model.pt model.safetensors
    semibin_model_convert.py safetensors2pt model.safetensors model.pt
"""
import json
import os
import sys
import tempfile

import torch
from safetensors import safe_open
from safetensors.torch import load_file, save_file


def to_safetensors(source, target):
    checkpoint = torch.load(source, map_location="cpu", weights_only=True)
    state_dict = checkpoint["model_state_dict"]
    save_file(
        state_dict,
        target,
        metadata={
            "model_name": checkpoint["model_name"],
            "params": json.dumps(checkpoint["params"]),
        },
    )
    check_round_trip(target, state_dict)


def to_pytorch(source, target):
    try:
        with safe_open(source, framework="pt", device="cpu") as handle:
            metadata = dict(handle.metadata() or {})
        state_dict = load_file(source, device="cpu")
    except Exception:
        # not a safetensors file: assume it is already a PyTorch checkpoint
        torch.save(torch.load(source, map_location="cpu", weights_only=True), target)
        return
    missing = {"model_name", "params"} - set(metadata)
    if missing:
        raise SystemExit(
            f"'{source}' lacks the {', '.join(sorted(missing))} metadata SemiBin needs to load the model"
        )
    torch.save(
        {
            "model_name": metadata["model_name"],
            "model_state_dict": state_dict,
            "params": json.loads(metadata["params"]),
        },
        target,
    )


def check_round_trip(path, state_dict):
    """Make sure that the model `SemiBin2 bin` reads back is the trained one.

    SemiBin needs the model name and the input dimension (both are stored as
    safetensors metadata) to rebuild the network, so a conversion that only
    keeps the weights would produce a file that cannot be used for binning.
    """
    from SemiBin.semi_supervised_model import model_load

    with safe_open(path, framework="pt", device="cpu") as handle:
        metadata = handle.metadata() or {}
    for key in ("model_name", "params"):
        if key not in metadata:
            raise SystemExit(f"'{path}' lacks the {key!r} metadata")
    with tempfile.TemporaryDirectory() as tmpdir:
        roundtrip = os.path.join(tmpdir, "model.pt")
        to_pytorch(path, roundtrip)
        model = model_load(roundtrip, "cpu")
    if set(model.state_dict()) != set(state_dict):
        raise SystemExit(f"'{path}' does not round trip to the trained weights")


def main():
    direction, source, target = sys.argv[1:4]
    if direction == "pt2safetensors":
        to_safetensors(source, target)
    elif direction == "safetensors2pt":
        to_pytorch(source, target)
    else:
        raise SystemExit("Unknown direction: " + direction)


if __name__ == "__main__":
    main()
