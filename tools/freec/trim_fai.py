chr_names = []
with open("capture.bed") as capture:
    for line in capture.readlines():
        chr_name = line.split()[0]
        if chr_name not in chr_names:
            chr_names.append(chr_name)

with open("genome.fa.fai") as fai:
    for line in fai.readlines():
        chr_name = line.split()[0]
        if chr_name in chr_names:
            with open("genome_trimmed.fa.fai", "a") as fai_trimmed:
                fai_trimmed.write(line)
