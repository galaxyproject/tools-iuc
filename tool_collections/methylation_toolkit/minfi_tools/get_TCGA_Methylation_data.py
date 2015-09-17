import sys
import tarfile


# Un tar file
if (tarfile.is_tarfile(sys.argv[1])):
    # print "File %s is a tarfile" %  sys.argv[1]
    tar = tarfile.open(sys.argv[1])
    tar.extractall()
    tar.close()
