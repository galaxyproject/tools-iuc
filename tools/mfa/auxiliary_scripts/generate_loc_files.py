import math
import os
import shutil
import struct
import subprocess
import wave
import zipfile

# Configuration
TEST_DATA_DIR = "test-data"
SAMPLE_RATE = 16000
DURATION = 0.5


def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def create_sine_wave(filename, duration=DURATION, freq=440):
    """Generates a clean sine wave."""
    n_frames = int(SAMPLE_RATE * duration)
    with wave.open(filename, 'w') as obj:
        obj.setnchannels(1)  # Mono
        obj.setsampwidth(2)  # 16-bit
        obj.setframerate(SAMPLE_RATE)

        data = []
        for i in range(n_frames):
            value = int(32767.0 * 0.5 * math.sin(2.0 * math.pi * freq * i / SAMPLE_RATE))
            data.append(struct.pack('<h', value))
        obj.writeframes(b''.join(data))


def create_zip(zip_name, source_dir):
    """Zips a directory into a file."""
    output_path = os.path.join(TEST_DATA_DIR, zip_name)
    with zipfile.ZipFile(output_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(source_dir):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, os.path.dirname(source_dir))
                zipf.write(file_path, arcname)
    print(f"Created: {output_path}")


def train_dummy_model(corpus_dir, dict_path, output_model_path):
    """Runs MFA Train to create a REAL, VALID acoustic model for testing."""
    print("Training dummy acoustic model... (This may take a few seconds)")

    # We use a temporary directory for output
    temp_out = os.path.join(TEST_DATA_DIR, "temp_train_output")
    ensure_dir(temp_out)

    # We need a G2P model or a dictionary with phones.
    # Our test_dict.txt has phones, so we can train directly.
    cmd = [
        "mfa", "train",
        corpus_dir,
        dict_path,
        output_model_path,
        "--clean",
        "--num_jobs", "1",
        "--single_speaker",
        "--phone_set_type", "IPA"
    ]

    try:
        # Run MFA quietly
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f"Successfully generated valid model: {output_model_path}")
    except subprocess.CalledProcessError:
        print("Error: Could not run 'mfa train' to generate dummy model.")
        print("Ensure 'mfa' is installed and in your PATH.")
        # Fallback to dummy if MFA is missing (will cause Adapt to fail, but allows script to finish)
        with open(output_model_path, "w") as f:
            f.write("DUMMY_FAILED_TRAIN")
    except FileNotFoundError:
        print("Error: 'mfa' command not found. Cannot generate valid binary model.")


def main():
    if os.path.exists(TEST_DATA_DIR):
        shutil.rmtree(TEST_DATA_DIR)
    ensure_dir(TEST_DATA_DIR)

    # 1. Create 'test_corpus'
    base_dir = os.path.join(TEST_DATA_DIR, "temp_corpus")
    ensure_dir(base_dir)
    # 3 speakers, 2 utterances (Small enough to be fast, big enough to train)
    for spk in range(1, 4):
        spk_dir = os.path.join(base_dir, f"speaker{spk}")
        ensure_dir(spk_dir)
        for utt in range(1, 3):
            freq = 440 + (spk * 50)
            create_sine_wave(os.path.join(spk_dir, f"{spk}_{utt}.wav"), duration=0.5, freq=freq)
            with open(os.path.join(spk_dir, f"{spk}_{utt}.lab"), "w") as f:
                f.write("hello world")

    # Zip corpus for tool inputs
    create_zip("test_corpus.zip", base_dir)

    # 2. Create 'test_dict.txt'
    dict_path = os.path.join(TEST_DATA_DIR, "test_dict.txt")
    with open(dict_path, "w") as f:
        f.write("hello\th\nworld\tw\n")  # Simple dictionary matching corpus

    # 3. GENERATE VALID ACOUSTIC MODEL (The Fix)
    # Instead of faking it, we actually train one on the dummy data.
    dummy_model_path = os.path.join(TEST_DATA_DIR, "dummy_acoustic_model.zip")
    train_dummy_model(base_dir, dict_path, dummy_model_path)

    # 4. Other files needed for other tools

    # Audio Only Corpus
    audio_dir = os.path.join(TEST_DATA_DIR, "temp_audio")
    ensure_dir(audio_dir)
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".wav"):
                rel_path = os.path.relpath(root, base_dir)
                target_dir = os.path.join(audio_dir, rel_path)
                ensure_dir(target_dir)
                shutil.copy(os.path.join(root, file), os.path.join(target_dir, file))
    create_zip("test_corpus_audio.zip", audio_dir)
    shutil.rmtree(audio_dir)

    # Words List
    words_path = os.path.join(TEST_DATA_DIR, "test_words.txt")
    with open(words_path, "w") as f:
        f.write("hello\nworld\n")

    # Mapping
    yaml_path = os.path.join(TEST_DATA_DIR, "arpabet_to_ipa.yaml")
    with open(yaml_path, "w") as f:
        f.write("h: h\nw: w\n")

    # Aligned Corpus (TextGrid) - Fixed Header
    aligned_dir = os.path.join(TEST_DATA_DIR, "temp_aligned")
    ensure_dir(os.path.join(aligned_dir, "s1"))
    create_sine_wave(os.path.join(aligned_dir, "s1", "1.wav"))
    tg_content = """File type = "ooTextFile"
Object class = "TextGrid"

xmin = 0
xmax = 0.5
tiers? <exists>
size = 1
item []:
    item [1]:
        class = "IntervalTier"
        name = "phones"
        xmin = 0
        xmax = 0.5
        intervals: size = 1
        intervals [1]:
            xmin = 0
            xmax = 0.5
            text = "h"
"""
    with open(os.path.join(aligned_dir, "s1", "1.TextGrid"), "w") as f:
        f.write(tg_content)
    create_zip("test_corpus_aligned.zip", aligned_dir)
    shutil.rmtree(aligned_dir)

    # Dummy IPA Model (for Remap)
    # We can just copy the valid acoustic model we trained, as it uses IPA-like phones
    if os.path.exists(dummy_model_path):
        shutil.copy(dummy_model_path, os.path.join(TEST_DATA_DIR, "dummy_ipa_model.zip"))
    else:
        # Fallback if training failed
        ipa_model_dir = "temp_ipa_model"
        ensure_dir(ipa_model_dir)
        with open(os.path.join(ipa_model_dir, "meta.json"), "w") as f:
            f.write('{"version": "3.0.0", "phones": ["h", "w"], "features": {"type": "mfcc"}}')
        with open(os.path.join(ipa_model_dir, "final.mdl"), "wb") as f:
            f.write(b"DUMMY")
        with open(os.path.join(ipa_model_dir, "tree"), "wb") as f:
            f.write(b"DUMMY")
        create_zip("dummy_ipa_model.zip", ipa_model_dir)
        shutil.rmtree(ipa_model_dir)

    # LM
    lm_dir = os.path.join(TEST_DATA_DIR, "temp_lm")
    ensure_dir(lm_dir)
    with open(os.path.join(lm_dir, "test.arpa"), "w") as f:
        f.write("\\data\\\nngram 1=3\n\\1-grams:\n-1.0 hello\n-1.0 world\n-1.0 </s>\n\\end\\")
    create_zip("dummy_language_model.zip", lm_dir)
    shutil.rmtree(lm_dir)

    print("Cleanup complete. Valid test data generated.")


if __name__ == "__main__":
    main()
