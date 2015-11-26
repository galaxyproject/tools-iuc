import sys
import re

def update(inputname):
    tool = open (inputname,'r')
    out = open (inputname+'.new','w')

    for s in tool:
        if "requirement type=\"package\" version=\"1.27\"" in s: #mothur tool version + 1
            s = s.replace("1.27","1.33")
            #print s
            out.write(s,)
        elif "tool id" in s: #header version + 1
            version = re.findall(r'\d+',s)
            #finds the version. This assumes that the version is the only int, and takes the second param
            newversion = int(version[1])+1
            newversion = str(newversion)
            s = s.replace(version[1],newversion)
            #print out, s
            out.write(s,)
        else:
            #print out, s,
            out.write(s,)

