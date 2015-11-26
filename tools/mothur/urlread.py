#script that hopefully automates some xml editting by reading the page for you
#v0.0: I think I'll just make it search pages and output the new params - that way at least you can ignore no changes
#USAGE: give the location of the folder containing the tools as the first and only parameter
import urllib2
import StringIO
import sys
import re
from os import listdir
from os import chdir
import autoup
import sourceparse


def getrev (name):
#print "Enter the current version of the tool"
#print "Enter the desired version of the tool"
    print "==========================================================="
    filename = name;
    #cut the .xml out of the name so we can get the webpage
    name = name[:-4]
    print "Attemping ", name
    try:
        response = urllib2.urlopen('http://www.mothur.org/wiki/' + str(name))
    except:
        print "Probably is not a wiki page for this tool, or it caught a non-xml file"
        return
    html = response.read()
    dog="".join(html)
#iterates through the html by line
    foundrev = None
    first = False
    #variable used to determine whether this tool can be automatically updated
    #notably, we check for the keyword "added"
    auto = True
    for s in StringIO.StringIO(dog):
        if "id=\"Revisions\"" in s or "id=\"Revision\"" in s:
            #find the revision column of the HTML
            #the next few lines give the revision tags
            foundrev = True
            first = True
            continue
        if foundrev:
            if "</ul>" in s:
                foundrev = False
                break
            #check for updates between the sepcified verions
            if first:
                s = s[8:]
                first = False
            else:
                s = s[9:]
            if "Bug Fix" not in s:
                print s,
            if "added" in s or "Added" in s:
                auto=False
     #if auto is still true at this point, we can reasonably inference that this tool is a trivial upgrade
    if auto:
        print name," determined to be a trivial upgrade."
        autoup.update(filename)


dir = sys.argv[1]
#change to this directory for trivial updates
try:
    chdir(dir)
except:
    print "Failed to chdir"
for f in  listdir(dir):
    getrev(f)
