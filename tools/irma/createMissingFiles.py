import glob
import os
import subprocess

dirPrefix = "resultDir/"
expectedSegments = {"A_MP": 7, "A_NP": 5, "A_HA": 4, "A_PB1": 2,
                    "A_PB2": 1, "A_NA": 6, "A_PA": 3, "A_NS": 8}


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
    return [segment for segment in expectedSegments.keys()
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


def writeEmptyAmendedFasta(identifier):
    #  irma names these files like: resultDir/amended_consensus/resultDir_<segNr>.fa
    open(dirPrefix + "amended_consensus/resultDir_" + str(expectedSegments[identifier]) + ".fa", 'x').close()


def samtoolsSortAllBam():
    for segment in expectedSegments.keys():
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
        writeEmptyAmendedFasta(segment)
    samtoolsSortAllBam()
