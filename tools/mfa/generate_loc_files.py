import os

# CONFIGURATION
# Ensure this matches exactly where your models live on the server
MFA_MODELS_ROOT = "/data/db/mfa/pretrained_models"
OUTPUT_DIR = "tool-data"

MODEL_MAPPINGS = {
    "acoustic":       "mfa_acoustic.loc",
    "dictionary":     "mfa_dictionary.loc",
    "g2p":            "mfa_g2p.loc",
    "ivector":        "mfa_ivector.loc",
    "language_model": "mfa_language_model.loc",
    "tokenizer":      "mfa_tokenizer.loc"
}

def format_name(model_id):
    name = model_id.replace("_", " ").title()
    name = name.replace(" Mfa", " (MFA)").replace(" Arpa", " (ARPA)").replace(" Cv", " (CommonVoice)")
    return name

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    print(f"Scanning models in: {MFA_MODELS_ROOT}")

    for folder_name, loc_filename in MODEL_MAPPINGS.items():
        search_path = os.path.join(MFA_MODELS_ROOT, folder_name)
        output_path = os.path.join(OUTPUT_DIR, loc_filename)

        if not os.path.exists(search_path):
            print(f"[SKIP] {folder_name} (Not found)")
            continue

        print(f"[GENERATE] {loc_filename}")
        
        with open(output_path, "w") as f:
            # Standard Header
            f.write("#<value>\t<dbkey>\t<name>\t<path>\n")
            
            files = sorted(os.listdir(search_path))
            for filename in files:
                if filename.startswith(".") or filename.endswith(".json"): continue
                
                model_id = os.path.splitext(filename)[0]
                display_name = format_name(model_id)
                
                # CRITICAL: This constructs the absolute path
                # e.g. /data/db/mfa/pretrained_models/acoustic/english_mfa.zip
                abs_path = os.path.join(search_path, filename)
                
                # Write 4 columns: ID, dbkey(default), Name, Path
                f.write(f"{model_id}\tdefault\t{display_name}\t{abs_path}\n")

if __name__ == "__main__":
    main()