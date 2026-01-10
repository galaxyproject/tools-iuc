import math
import os
import random
import shutil
import struct
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
    """Generates a sine wave mixed with WHITE NOISE to prevent LDA/MLLT crashes."""
    n_frames = int(SAMPLE_RATE * duration)
    with wave.open(filename, 'w') as obj:
        obj.setnchannels(1)  # Mono
        obj.setsampwidth(2)  # 16-bit
        obj.setframerate(SAMPLE_RATE)

        data = []
        for i in range(n_frames):
            tone = 1000 * math.sin(2.0 * math.pi * freq * i / SAMPLE_RATE)
            noise = random.randint(-200, 200)
            value = int(tone + noise)
            value = max(-32767, min(32767, value))
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


def main():
    if os.path.exists(TEST_DATA_DIR):
        shutil.rmtree(TEST_DATA_DIR)
    ensure_dir(TEST_DATA_DIR)

    # 1. Create 'test_corpus'
    base_dir = "temp_corpus"
    ensure_dir(base_dir)
    for spk in range(1, 6):
        spk_dir = os.path.join(base_dir, f"speaker{spk}")
        ensure_dir(spk_dir)
        for utt in range(1, 4):
            freq = 440 + (spk * 50)
            create_sine_wave(os.path.join(spk_dir, f"{spk}_{utt}.wav"), duration=1.0, freq=freq)
            with open(os.path.join(spk_dir, f"{spk}_{utt}.lab"), "w") as f:
                f.write("hello world")
    create_zip("test_corpus.zip", base_dir)

    # 2. Create 'test_corpus_audio'
    audio_dir = "temp_audio"
    ensure_dir(audio_dir)
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".wav"):
                rel_path = os.path.relpath(root, base_dir)
                target_dir = os.path.join(audio_dir, rel_path)
                ensure_dir(target_dir)
                shutil.copy(os.path.join(root, file), os.path.join(target_dir, file))
    create_zip("test_corpus_audio.zip", audio_dir)

    # 3. Create 'test_corpus_raw'
    raw_dir = "temp_raw"
    ensure_dir(os.path.join(raw_dir, "jp_speaker"))
    create_sine_wave(os.path.join(raw_dir, "jp_speaker", "jp1.wav"), freq=500)
    with open(os.path.join(raw_dir, "jp_speaker", "jp1.lab"), "w") as f:
        f.write("こんにちは世界")
    create_zip("test_corpus_raw.zip", raw_dir)

    # 4. Create 'test_dict.txt'
    dict_path = os.path.join(TEST_DATA_DIR, "test_dict.txt")
    with open(dict_path, "w") as f:
        f.write("hello\tHH AH0 L OW1\n")
        f.write("world\tW ER1 L D\n")
        f.write("yes\tY EH1 S\n")
        f.write("no\tN OW1\n")

    # 5. Create 'test_words.txt'
    words_path = os.path.join(TEST_DATA_DIR, "test_words.txt")
    with open(words_path, "w") as f:
        f.write("test\n")

    # 6. Create 'arpabet_to_ipa.yaml'
    yaml_path = os.path.join(TEST_DATA_DIR, "arpabet_to_ipa.yaml")
    with open(yaml_path, "w") as f:
        f.write("HH: h\n")
        f.write("AH0: ə\n")
        f.write("L: l\n")
        f.write("OW1: oʊ\n")
        f.write("W: w\n")
        f.write("ER1: ɝ\n")
        f.write("D: d\n")
        f.write("Y: j\n")
        f.write("EH1: ɛ\n")
        f.write("S: s\n")
        f.write("N: n\n")

    # 7. Create 'dummy_acoustic_model.zip' (ARPABET)
    model_dir = "temp_model"
    ensure_dir(model_dir)
    with open(os.path.join(model_dir, "meta.json"), "w") as f:
        f.write('{"version": "3.0.0", "phones": ["HH", "AH0", "L", "OW1", "W", "ER1", "D", "Y", "EH1", "S", "N"], "features": {"type": "mfcc"}}')
    with open(os.path.join(model_dir, "final.mdl"), "wb") as f:
        f.write(b"DUMMY")
    with open(os.path.join(model_dir, "tree"), "wb") as f:
        f.write(b"DUMMY")
    create_zip("dummy_acoustic_model.zip", model_dir)
    shutil.rmtree(model_dir)

    # 8. Create 'test_corpus_aligned.zip'
    aligned_dir = "temp_aligned"
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
            text = "HH"
"""
    with open(os.path.join(aligned_dir, "s1", "1.TextGrid"), "w") as f:
        f.write(tg_content)
    create_zip("test_corpus_aligned.zip", aligned_dir)
    shutil.rmtree(aligned_dir)

    # 9. Create 'dummy_ipa_model.zip' (TARGET PHONES for Remap Test)
    # This model defines the IPA phones we are mapping TO.
    ipa_model_dir = "temp_ipa_model"
    ensure_dir(ipa_model_dir)
    with open(os.path.join(ipa_model_dir, "meta.json"), "w") as f:
        f.write('{"version": "3.0.0", "phones": ["h", "ə", "l", "oʊ", "w", "ɝ", "d", "j", "ɛ", "s", "n"], "features": {"type": "mfcc"}}')
    with open(os.path.join(ipa_model_dir, "final.mdl"), "wb") as f:
        f.write(b"DUMMY")
    with open(os.path.join(ipa_model_dir, "tree"), "wb") as f:
        f.write(b"DUMMY")
    create_zip("dummy_ipa_model.zip", ipa_model_dir)
    shutil.rmtree(ipa_model_dir)

    # 10. Create 'dummy_language_model.zip' (For Transcribe)
    lm_dir = "temp_lm"
    ensure_dir(lm_dir)
    with open(os.path.join(lm_dir, "test.arpa"), "w") as f:
        f.write("\\data\\\nngram 1=3\n\\1-grams:\n-1.0 hello\n-1.0 world\n-1.0 </s>\n\\end\\")
    create_zip("dummy_language_model.zip", lm_dir)
    shutil.rmtree(lm_dir)

    # Cleanup
    shutil.rmtree(base_dir)
    shutil.rmtree(audio_dir)
    shutil.rmtree(raw_dir)
    print("Cleanup complete. Test data generated.")


if __name__ == "__main__":
    main()
