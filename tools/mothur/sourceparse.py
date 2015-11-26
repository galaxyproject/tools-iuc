import sys
from itertools import izip
import re

#get file A (old date)
#get file B (new ver)
#grep both for "CommandParam"
#sort both files
#compare line by line, flag "new" lines
#ouput only these new lines

def getSourceDiff (filea,fileb):

#get the file names
    old = filea
    new = fileb

#open these things
    oldread = open(old,'r')
    newread = open(new,'r')
#new files
    oldout = open(old+'.new','w+')
    newout= open(new+'.new','w+')


#grep&sort both files, write to new
    def grepSort(fin,fout):
        mylist = []
        #im guessing these things are strings or something, grep these?
        for l in fin:
            if "CommandParameter" in l:
                mylist.append(l)
            

        mylist.sort()

        for item in mylist:
            print>>fout, item,

    grepSort(oldread,oldout)
    grepSort(newread,newout)

    oldread.close()
    newread.close()

#compare the files line by line. 

    newout.seek(0,0)
    oldout.seek(0,0)

    while True:
        line1 = oldout.readline()
        line2 = newout.readline()

        if line1 == "" and line2 != "":
            #reached end of old file, get new stuff from line2
            print line2,

        if not line2:
            break


