import glob
import os
import subprocess

dirPrefix = "resultDir/"
expectedSegments = ["A_MP", "A_NP", "A_HA", "A_PB1",
                    "A_PB2", "A_NA", "A_PA", "A_NS"]


def renameSubtypeFiles(identifier):
    files = glob.glob(dirPrefix + "A_" + identifier + "_*.*")
    for file in files:
        ext = file.split('.')[-1]
        os.rename(file, dirPrefix + "A_" + identifier + "." + ext)


def getMissingSegments():
    presentSegments = []
    for file in os.listdir(dirPrefix):
        if file.endswith(".fasta"):
            presentSegments.append(file.split('.')[0])
    return [segment for segment in expectedSegments
            if segment not in presentSegments]


def getBamHeaderFromAnyFile():
    anyBamFile = glob.glob(dirPrefix + "*.bam")[0]
    samtoolsCmd = ["samtools", "view", "-H", anyBamFile]
    result = subprocess.check_output(samtoolsCmd, text=True)
    return result.split('\n')[0]


def getVcfHeaderFromAnyFile():
    with open(glob.glob(dirPrefix + "*.vcf")[0]) as f:
        anyVersionAndDateLines = f.readline() + f.readline()
        emptyHeaderLine = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        return anyVersionAndDateLines + emptyHeaderLine


def writeEmptyBam(identifier, bamHeader):
    with open("headerSamFile.sam", "w") as f:
        f.write(bamHeader)  # write header to a temporary sam file
    cmd = ['samtools', 'view', '-H', '-b', 'headerSamFile.sam']
    targetBam = dirPrefix + identifier + ".bam"
    with open(targetBam, "xb") as tB:
        subprocess.check_call(cmd, stdout=tB)
        os.remove("headerSamFile.sam")


def writeEmptyFasta(identifier):
    open(dirPrefix + identifier + ".fasta", 'x').close()


def writeEmptyVcf(identifier, vcfHeader):
    with open(dirPrefix + identifier + ".vcf", 'x') as f:
        f.write(vcfHeader)


def samtoolsSortAllBam():
    for segment in expectedSegments:
        os.rename(dirPrefix + segment + ".bam",
                  dirPrefix + segment + "_unsorted.bam")
        cmd = ['samtools', 'sort', dirPrefix + segment + "_unsorted.bam"]
        targetBam = dirPrefix + segment + ".bam"
        with open(targetBam, "w") as tB:
            subprocess.check_call(cmd, stdout=tB, text=True)
            os.remove(dirPrefix + segment + "_unsorted.bam")


if __name__ == "__main__":
    renameSubtypeFiles("HA")
    renameSubtypeFiles("NA")
    bamHeader = getBamHeaderFromAnyFile()
    vcfHeader = getVcfHeaderFromAnyFile()
    for segment in getMissingSegments():
        writeEmptyBam(segment, bamHeader)
        writeEmptyFasta(segment)
        writeEmptyVcf(segment, vcfHeader)
    samtoolsSortAllBam()
