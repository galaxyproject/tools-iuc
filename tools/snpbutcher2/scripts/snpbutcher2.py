#!/usr/bin/env python2

import sys, os.path
from commands import getoutput

#Clear logfile
logfile=open('autoconfig_snpbutcher.txt','w')


name = sys.argv[0]

def usage():
    print >> sys.stderr, '''
Usage: %s <map file> <maf file> <genotypes file> [OPTIONS]

   emmmmmm~~~~~~~~~~oT                             ver 2.8
   `""""""|          |
          |          |
          `----------'          
 OPTIONS:
  --hg19,     use human genome ver 19 (instead of the version 18)

  -u or --uninformative,     to include uninformative markers too
  -z or --zeromaf, to include zero minor allele frequency markers

  -nh or --nohapmap, to not check hapmap data

  -cm N or --centimorgans N  space markers at N cM (default=0.05)

  -s N or --spacing N   OR   -m N or --max N  
      select markers at intervals of N, or specify the max number
      of markers N, or select at intervals of N centiMorgans
      Both -u and -z are autoselected for this flag.
    ''' % name

    sys.exit(-1)

#===================
def parse_args():
     lenarg = len(sys.argv)

     hg19_flag = False
     uninform = False
     zeromaf=False
     usehapmap=True
     interval = max_interval = centim = -1

     files = []
     index = lenarg-1

     while index > 0:
        arg = sys.argv[index]

        if arg[0]=='-':

           if arg=='--hg19':
               hg19_flag=True
               sys.argv.remove(arg)

           elif arg=='-nh' or arg=="--nohapmap":
               usehapmap=False
               sys.argv.remove(arg)

           elif arg=='-u' or arg=="--uninformative":
               uninform=True
               sys.argv.remove(arg)

           elif arg=='-z' or arg=="--zeromaf":
               zeromaf=True
               sys.argv.remove(arg)

           elif arg=='-s' or arg=='--spacing':
               try:
                  interval=sys.argv[index+1]
                  sys.argv.remove(arg)
                  sys.argv.remove(interval)
                  interval=int(interval)

                  if interval < 0:
                     print >> sys.stderr,"[Error] Invalid Interval"
                     exit()
               except ValueError:
                  print >> sys.stderr, "[Error] Spacing: Not a number"
                  exit()
               except IndexError:
                  print >> sys.stderr,"[Error] Must specify an interval amount!"     
                  exit()

           elif arg=='-m' or arg=='--max':
               try:
                  max_interval=sys.argv[index+1]
                  sys.argv.remove(arg)
                  sys.argv.remove(max_interval)
                  max_interval=int(max_interval)

                  if max_interval < 0:
                     print >> sys.stderr,"[Error] Invalid Max Interval"
                     exit()
               except ValueError:
                  print >> sys.stderr, "[Error] Max: Not a number"
                  exit()
               except IndexError:
                  print >> sys.stderr,"[Error] Must specify a max amount!"     
                  exit()

           elif arg=='-cm' or arg=='--centim' or '--centimorgans':
               try:
                  centim=sys.argv[index+1]
                  sys.argv.remove(arg)
                  sys.argv.remove(centim)
                  centim=float(centim)

                  if centim < 0:
                     print >> sys.stderr,"[Error] Invalid CentiMorgan Interval"
                     exit()
               except ValueError:
                  print >> sys.stderr, "[Error] CentiMorgan: Not a number"
                  exit()
               except IndexError:
                  print >> sys.stderr,"[Error] Must specify a centimorgan amount!"     
                  exit()

           else:
               usage()

        index -= 1

     
     sys.argv.remove(name)
     files=sys.argv

     if len(files)<3:
         usage()

     if ((interval!=-1 and max_interval!=-1)
	 or (interval!=-1 and centim!=-1) 
	 or (centim!=-1 and max_interval!=-1) ):
         print >> sys.stderr, "[Error] Cannot specify more than one interval type"
         sys.exit(-1)


     # Spacing arguments assume zeromaf and use uninf
     if interval!=-1 or max_interval!=-1:
          zeromaf = uninform = True

     print >> sys.stderr, "[ARGS]",

     if max_interval!=-1:print >> sys.stderr, "Max: %d," % max_interval,
     if interval!=-1:print >> sys.stderr, "Spacing: %d," % interval,
     if centim!=-1:print >> sys.stderr, "cM: %f," % centim,
     if hg19_flag:print >> sys.stderr, "Human Genome v19,",
     if uninform:print >> sys.stderr, "Uninformative too,",
     if zeromaf:print >> sys.stderr, "Zero maf markers too,",
     if not(usehapmap):print >> sys.stderr, "NOT using hapmap",
  
     print >> sys.stderr, "\n"

     print >> logfile,\
"\ninterval=",interval,\
"\nmax_interval=",max_interval,\
"\nusehapmap=",usehapmap,\
"\nhg19_flag=",hg19_flag,\
"\ncM=",centim,\
"\nzeromaf=",zeromaf,\
"\nuninformative=",uninform
     return files, interval, max_interval, hg19_flag, uninform, zeromaf, centim, usehapmap

#=====================================================

def uninformativeSNPs(genotypesfile):
	uninformative = {};
	genotypecount = 0


	print >> sys.stderr, "[INFO] reading '%s'..." % genotypesfile

	f = open(genotypesfile)
	for line in f :
#		if line.startswith('rs') :
		    tmp = line.strip().split()
		    marker = tmp[0]
		    gtypes = tmp[1:]
		    gtypes = set(filter(lambda x : x != 'NC', gtypes))
		    uninformative[marker] = True if len(gtypes) == 1 else False
#		genotypecount += 1
		    genotypecount += 1
	f.close()

	print >> sys.stderr, "[INFO] found %d/%d SNPs uninformative" % (len(filter(lambda x : x == True, uninformative.values())), genotypecount)
	
	return uninformative

#==================

def minorAlleleFreqs(maffile):
	frequencies = {}

	print >> sys.stderr, "[INFO] reading '%s'..." % maffile

	f = open(maffile)
	f.readline() #strip header
	for line in f :
#		if line.startswith('rs') :
		if len(line)<5:continue
#		try:
		rs,freq = line.strip().split()
		frequencies[rs] = float(freq)
#		except ValueError:
#			print >> sys.stderr, '\n', line
#			exit(-1)
	f.close()

	print >> sys.stderr, "[INFO] read %d SNP minor allele frequencies" % len(frequencies)

	return frequencies

#==================

def get_hapmap_map(i) :
    data = {}
    #print >> sys.stderr, "[INFO] reading c%d reference..." % i
    global hg19_flag

    try:
       i = int(i)
    except ValueError:
       pass

    i = str(i)
        
    try:
       if hg19_flag:
          g = open(hapmap_path + "genetic_map_GRCh37_chr%s.txt" % i)
       else:
          g = open(hapmap_path + "genetic_map_chr%s_b36.txt" % i)
    except IOError:
#       error_msg = "No Hapmap data for chr"+i+(" hg19" if hg19_flag else " hg18")+" - using uncorrected cM positions in original map"
 #      print >> sys.stderr, "\n[Caution]",error_msg,
  #     print >> logfile, "error: "+error_msg		#log message in readme_inp

       return -1
    
    g.readline()
    
    if hg19_flag:
        for line in g :
           chr_nm,hm_phy,hm_com,hm_gen = line.strip().split()
           data[hm_phy] = float(hm_gen)

    else:
        for line in g :
           hm_phy,hm_com,hm_gen = line.strip().split()
           data[hm_phy] = float(hm_gen)
    
    g.close()
    
    return data

#==================

# Default method of finding all the informative SNPs in the map file
def informativeSNPs(mapfile, uninformative, frequencies, 
	spacing, zeromaf, use_uninf, checkhap, interval=0.05):

	geneticmap_data = None
	markercount = 0
	threshold=0.0

	print >> sys.stderr, "[INFO] reading '%s' and printing new map..." % mapfile
	print >> sys.stderr, "[INFO] chromosome",

	f = open(mapfile)

	print f.readline().strip()
	current_chr = 0
	line_count=0;

	for line in f :
		chr,rsid,genpos,physpos,rsid2,junk = line.strip().split()
		genpos = float(genpos)
#		try :
#			chr = int(chr)
#		except :
			#pass # X, XY, Y
#			print line.strip()
#			continue
		
		if chr != current_chr :
			current_chr = chr
			threshold = 0.0
			print >> sys.stderr, "%s" % chr,
			geneticmap_data = get_hapmap_map(chr)

		skip=False

		if checkhap:
			if geneticmap_data != -1: 		#Has hapmap dataa:
				if geneticmap_data.has_key(physpos) :
					genpos = geneticmap_data[physpos]	#Update genpos with hapmap data
				else:
					skip = True				#Has hapmap, no SNP? Skip


		if not(skip):
			if spacing!=-1:
				if line_count%spacing==0:
					print '\t'.join([chr, rsid, str(genpos), physpos, rsid, junk])
					markercount += 1
			else:		
				if genpos > threshold :
					try :
#						if not rsid.startswith("rs") :  # not kgp friendly this
#							continue
				    
						if not(use_uninf):
							if uninformative[rsid] :
								continue
				    
						if not(zeromaf):
							if (frequencies[rsid] == 0) or (frequencies[rsid] == 1) :
								continue

					except KeyError :
						continue
				
					threshold = genpos + interval

					print '\t'.join([chr, rsid, str(genpos), physpos, rsid, junk])
					markercount += 1		

		line_count += 1;

	f.close()
	print >> sys.stderr, "\n[INFO] new map contains %d SNPs" % markercount


###############################  M  A  I  N  ###################################

try:
	files, equally_spaced, max_interval, hg19_flag, use_uninf, zeromaf, cm, checkhap = parse_args()

	mapfile = files[0]
	maffile = files[1]
	genotypesfile = files[2]

	hapmap_path = '/usr/local/maps/hapmap_geneticmaps/'

	for i in [mapfile, maffile, genotypesfile] :
		if not os.path.exists(i) :
		    print >> sys.stderr, "[ERROR] '%s' does not exist" % i
		    sys.exit(-1)

	frequencies = [] if zeromaf else minorAlleleFreqs(maffile)
	
	uninformative = [] if use_uninf else uninformativeSNPs(genotypesfile)
	#print >> sys.stderr, uninformative

	# Parse Spacing args
	spacing = -1

	if equally_spaced!=-1:
		spacing = equally_spaced
		spacing = 1 if (spacing < 1) else spacing

	if max_interval!=-1:
		num_lines_mapf = int(getoutput("wc -l "+mapfile+" | awk '{print $1}'"))
		spacing = int( float(num_lines_mapf) / float(max_interval) )-1
		spacing = 1 if (spacing < 1) else spacing
	

	# Parse centimorgan spacing
	interval=0.05 if (cm < 0) else cm
	mark = informativeSNPs(mapfile, uninformative, frequencies, spacing, zeromaf, use_uninf, checkhap, interval )
	
################################ L O G G I N G #################################
### Added config output for readme_gen.py ###  mct 2013/01/11
	#Log error for readme_gen.py to print in README.txt
	print >> logfile, "ver:", ("hg19" if hg19_flag else "hg18")
################################################################################
except KeyboardInterrupt:
   print >> sys.stderr, "\r[Terminated]"
   exit(-1)
