import sys
import tarfile


# Un tar file
if (tarfile.is_tarfile(sys.argv[1])):
    # print "File %s is a tarfile" %  sys.argv[1]
    tar = tarfile.open(sys.argv[1])
    tar.extractall()
    tar.close()

# print str(sys.argv[1])
# folder_name = (str(sys.argv[1])).replace(".tar.gz", "")
# print folder_name
# os.rename(folder_name, "gdac_methylation_dir")
# for filename in os.listdir("gdac_methylation_dir"):
#     if os.path.getsize(filename) > 3000:
#         os.rename(filename, "METHYLATION_DATA.txt")
