#!/bin/bash
# Download CoreNLP English models for testing
MODELS_VERSION="4.5.10"
MODELS_URL="https://repo1.maven.org/maven2/edu/stanford/nlp/stanford-corenlp/${MODELS_VERSION}/stanford-corenlp-${MODELS_VERSION}-models-english.jar"
MODELS_FILE="stanford-corenlp-${MODELS_VERSION}-models-english.jar"

if [ ! -f "${MODELS_FILE}" ]; then
    echo "Downloading CoreNLP English models..."
    curl -L -o "${MODELS_FILE}" "${MODELS_URL}"
fi
