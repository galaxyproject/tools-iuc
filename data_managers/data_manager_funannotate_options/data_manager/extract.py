import json
import os
import sys

fun_db = sys.argv[1]
fun_db_value = sys.argv[2]
dmjson = sys.argv[3]

content = []
# get options for parameter --busco_db
# which are just the subfolders of the db dir (minus outgroups/ and trained_species/)
# https://github.com/nextgenusfs/funannotate/blob/8cc40728fee61566fdf736c1f2292e14cc117660/funannotate/predict.py#L319
for d in os.scandir(fun_db):
    if not d.is_dir():
        continue
    if d.name in ['outgroups', 'trained_species']:
        continue
    if not os.path.exists(os.path.join(d, "dataset.cfg")):
        continue
    name = d.name.replace("_", " ").capitalize()
    content.append({'value': d.name, 'name': name, 'select': 'busco_db', 'db_value': fun_db_value})

# --busco_seed_species
# trained_species
for d in os.scandir(os.path.join(fun_db, "trained_species")):
    if not d.is_dir():
        continue
    if not os.path.exists(os.path.join(d, "info.json")):
        continue
    name = d.name.replace("_", " ").capitalize()
    content.append({'value': d.name, 'name': name, 'select': 'trained_species', 'db_value': fun_db_value})

# --busco_seed_species
# outgroups
for f in os.scandir(os.path.join(fun_db, "outgroups")):
    if f.is_dir():
        continue
    if not f.name.endswith("_buscos.fa"):
        continue
    value = f.name[:-10]
    name = ' - '.join([x.replace("_", " ").capitalize() for x in value.split('.')])
    content.append({'value': value, 'name': name, 'select': 'outgroup', 'db_value': fun_db_value})

with open(dmjson, "w") as fh:
    json.dump({"data_tables": {"funannotate_options": content}}, fh)

print(f'{len([c for c in content if c["select"]=="busco_db"])} x busco_db\n')
print(f'{len([c for c in content if c["select"]=="trained_species"])} x trained_species\n')
print(f'{len([c for c in content if c["select"]=="outgroup"])} x outgroup\n')
