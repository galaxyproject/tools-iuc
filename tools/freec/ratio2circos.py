import math
import sys

ploidy = int(sys.argv[1])

with open("./output/sample.bam_ratio.BedGraph") as bed:
    with open("./output/sample.bam_ratio_log2_circos.txt", "w+") as olog2r:
        for line in bed.readlines():
            ls = line.split()
            if ls[0] != "track" and float(ls[3]) > 0:
                log2_ratio = math.log2(float(ls[3]) / ploidy)
                olog2r.write("{}\t{}\t{}\t{}\n".format(ls[0], ls[1], ls[2], log2_ratio))

with open("./genome.fa.fai") as fai:
    with open("./output/karyotype_circos.txt", "w+") as ochr:
        for line in fai.readlines():
            ls = line.split()
            ochr.write("chr - {}\t{}\t0\t{}\t{}\n".format(ls[0], ls[0].strip("chr").lower(), ls[1], ls[0]))
