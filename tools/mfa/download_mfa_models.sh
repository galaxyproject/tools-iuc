#!/bin/bash

# =============================================================================
# Script to download all pretrained models for Montreal Forced Aligner (MFA)
# =============================================================================
#
# Usage:
#   bash download_mfa_models.sh [GITHUB_TOKEN]
#
# Arguments:
#   GITHUB_TOKEN: A GitHub personal access token to increase API 
#                            rate limits. Highly recommended for bulk downloads.
#
# =============================================================================

TARGET_DIR="/data/db/mfa"
export MFA_ROOT_DIR="$TARGET_DIR"

# 1. Check for Token
if [ -z "$GITHUB_TOKEN" ]; then
    if [ -n "$1" ]; then
        export GITHUB_TOKEN="$1"
        echo "Token detected from argument."
    else
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "WARNING: No GITHUB_TOKEN detected."
        echo "You will likely hit API rate limits immediately."
        echo "Usage: export GITHUB_TOKEN='...' && bash script.sh"
        echo "Or:    bash script.sh 'GITHUB_TOKEN'"
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "Waiting 5 seconds before proceeding (Ctrl+C to cancel)..."
        sleep 5
    fi
else
    echo "GitHub Token detected in environment."
fi

# 2. Check Permissions
if [ ! -w "$TARGET_DIR" ]; then
    echo "Error: Directory $TARGET_DIR is not writable."
    exit 1
fi

echo "==================================================="
echo "MFA Downloader V4"
echo "Storage Location: $MFA_ROOT_DIR"
echo "==================================================="

MODEL_TYPES=("acoustic" "dictionary" "g2p" "ivector" "language_model" "tokenizer")

for type in "${MODEL_TYPES[@]}"; do
    echo "---------------------------------------------------"
    echo "Fetching list for model type: '$type'..."
    
    # Capture output
    RAW_OUTPUT=$(mfa model download "$type" 2>&1)
    
    # Parse keys using Python
    MODELS=$(echo "$RAW_OUTPUT" | python3 -c "import sys, re; txt = sys.stdin.read(); print('\n'.join(re.findall(r\"'([a-zA-Z0-9_\-\.]+)':\", txt)))")
    
    # ERROR HANDLING: If no models found, show the raw output to diagnose WHY
    if [ -z "$MODELS" ]; then
        echo "  [WARNING] Could not parse model list for '$type'."
        echo "  --- RAW OUTPUT FROM MFA ---"
        echo "$RAW_OUTPUT"
        echo "  ---------------------------"
        continue
    fi
    
    # Download Loop
    while IFS= read -r model_name; do
        if [ -n "$model_name" ]; then
            echo -n "  -> Downloading: $model_name ... "
            
            # Run download, suppress output but capture exit code
            mfa model download "$type" "$model_name" --ignore_cache > /dev/null 2>&1
            
            if [ $? -eq 0 ]; then
                echo "[OK]"
            else
                echo "[FAIL]"
                # If a download fails, it might be rate limit.
                # Uncomment the line below to see why it failed:
                # mfa model download "$type" "$model_name" --ignore_cache
            fi
            
            # Sleep 1 second to be polite to the API
            sleep 1
        fi
    done <<< "$MODELS"
done

echo "==================================================="
echo "Process Finished."
