import math
import sys

ploidy = int(sys.argv[1])

with open("./output/sample.bam_ratio.BedGraph") as bed:
    for line in bed.readlines():
        ls = line.split()
        if ls[0] != "track" and float(ls[3]) > 0:
            log2_ratio = math.log2(float(ls[3]) / ploidy)
            with open("./output/sample.bam_ratio_log2_circos.txt", "a") as olog2r:
                olog2r.write("{}\t{}\t{}\t{}\n".format(ls[0], ls[1], ls[2], log2_ratio))

with open("./genome.fa.fai") as fai:
    for line in fai.readlines():
        ls = line.split()
        with open("./output/karyotype_circos.txt", "a") as ochr:
            ochr.write("chr - {}\t{}\t0\t{}\t{}\n".format(ls[0], ls[0].strip("chr").lower(), ls[1], ls[0]))
