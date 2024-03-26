inconf = open("config.json", "r").readlines()
with open("config.json.local", "w") as bak:
    bak.write("".join(inconf))
urlbase = "https://galaxy.genomicsvl-students.cloud.edu.au/jbrowse/hum/"
utag = '"uri":'
for i, row in enumerate(inconf):
    ispath = False
    if row.strip().startswith(utag):
        ispath = True
        parth = row.split(utag)[1].strip().replace('"', "").replace("'", "")
        inconf[i] = '%s "%s%s"' % (utag, urlbase, parth)
with open("config.json", "w") as outconf:
    outconf.write("".join(inconf))
