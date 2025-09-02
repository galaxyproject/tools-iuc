import os
import sys
import json

workdir = sys.argv[1]

for fname in os.listdir(workdir):
    if fname.endswith(".json"):
        filepath = os.path.join(workdir, fname)
        file_id = fname.split(".")[0]  

        with open(filepath, "r") as f:
            data = json.load(f)

        if data.get("id") == "tbprofiler":
            data["id"] = file_id
            with open(filepath, "w") as f:
                json.dump(data, f, indent=4)

            print(f"Updated {fname} -> id = {file_id}")
