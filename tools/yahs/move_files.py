import os
import shutil

files = os.listdir()
for file in files:
    if "yahs_out" in file and "final" in file:
        shutil.move(file, "final_outs/" + file)
    elif "yahs_out_init" in file:
        shutil.move(file, "initial_break/" + file)
    elif "_break.agp" in file:
        shutil.move(file, "agp_break/" + file)
    elif "yahs_out" in file and ".agp" in file:
        shutil.move(file, "agp_out/" + file)
