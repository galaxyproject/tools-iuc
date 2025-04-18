import json
import os
import re
import sys

busco_db = os.path.join(sys.argv[1], "lineages")
busco_db_value = sys.argv[2]
dmjson = sys.argv[3]

content = []
for d in os.scandir(busco_db):
    if not d.is_dir():
        continue
    if not os.path.exists(os.path.join(d, "dataset.cfg")):
        continue
    name = re.sub(r"_odb\d+", "", d.name)
    name = name.replace("_", " ").capitalize()
    content.append({'value': d.name, 'name': name, 'db_value': busco_db_value})

with open(dmjson, "w") as fh:
    json.dump({"data_tables": {"busco_database_options": content}}, fh)

print(f'{len(content)} x busco_db\n')
