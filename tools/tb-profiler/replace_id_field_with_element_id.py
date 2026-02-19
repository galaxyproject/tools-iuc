import json
import os
import re
import sys


input_file = sys.argv[1]
output_file = sys.argv[2]

safe_file_name = re.sub(r'[^\w\-_\.]', '_', os.path.basename(output_file))

if not safe_file_name.endswith(".results.json"):
    safe_file_name += ".results.json"

file_id = os.path.splitext(safe_file_name)[0]

with open(input_file, "r") as f:
    data = json.load(f)

if data.get("id") == "tbprofiler":
    data["id"] = file_id

with open(safe_file_name, "w") as f:
    json.dump(data, f, indent=4)

print(f"Copied {input_file} -> {safe_file_name}, updated id = {file_id}")
