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
        print("Rename:")
        print(file, dirPrefix + "A_" + identifier + "." + ext)


def getMissingSegments():
    presentSegments = [filename.split('.')[0].split('/')[-1]
                       for filename in glob.glob(dirPrefix + "*.fasta")]
    return [segment for segment in expectedSegments
            if segment not in presentSegments]


def getBamHeaderFromAnyFile():
    anyBamFile = glob.glob(dirPrefix + "*.bam")[0]
    samtoolsCmd = "samtools view -H " + anyBamFile
    result = subprocess.run(samtoolsCmd, shell=True,
                            stdout=subprocess.PIPE, text=True)
    print("Sample Header: (from: ", anyBamFile)
    print(result)
    header = result.stdout.split('\n')[0]
    return header


def writeEmptyBam(identifier, bamHeader):
    with open("headerSamFile.sam", "w") as f:
        f.write(bamHeader)  # write header to a temporary sam file
    cmd = "samtools view -H -b headerSamFile.sam > " \
          + dirPrefix + identifier + ".bam"  # convert to bam
    print("Create bam file: " + cmd)
    if subprocess.run(cmd, shell=True, text=True).returncode == 0:
        os.remove("headerSamFile.sam")  # delete temporary sam file
    else:
        raise ValueError("Could not create bam file!")


def writeEmptyFile(identifier, ext):
    open(dirPrefix + identifier + "." + ext, 'a').close()


def samtoolsSortAllBam():
    for segment in expectedSegments:
        os.rename(dirPrefix + segment + ".bam",
                  dirPrefix + segment + "_unsorted.bam")
        cmd = "samtools sort " + dirPrefix + segment + "_unsorted.bam > " \
              + dirPrefix + segment + ".bam"
        if subprocess.run(cmd, shell=True, text=True).returncode == 0:
            print("Sorted bam created. Original bam removed for: " + segment)
            os.remove(dirPrefix + segment + "_unsorted.bam")
        else:
            raise ValueError("Could not sort bam file!")


if __name__ == "__main__":
    renameSubtypeFiles("HA")
    renameSubtypeFiles("NA")
    for segment in getMissingSegments():
        writeEmptyBam(segment, getBamHeaderFromAnyFile())
        writeEmptyFile(segment, "fasta")
        writeEmptyFile(segment, "vcf")
    samtoolsSortAllBam()
