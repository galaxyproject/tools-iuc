#!/bin/bash

# Install spaCy model script
# Usage: install_model.sh <model_name>

set -e

MODEL_NAME="$1"

if [ -z "$MODEL_NAME" ]; then
    echo "Usage: $0 <model_name>"
    exit 1
fi

echo "Installing spaCy model: $MODEL_NAME"

# Check if model is already installed
if python -c "import spacy; spacy.load('$MODEL_NAME')" 2>/dev/null; then
    echo "Model $MODEL_NAME is already installed"
    exit 0
fi

# Try to download the model
echo "Downloading model $MODEL_NAME..."
python -m spacy download "$MODEL_NAME"

# Verify installation
if python -c "import spacy; spacy.load('$MODEL_NAME')" 2>/dev/null; then
    echo "Model $MODEL_NAME installed successfully"
else
    echo "Failed to install model $MODEL_NAME"
    exit 1
fi