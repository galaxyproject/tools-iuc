#!/usr/bin/env python

"""
http://www.mothur.org/

Supports mothur version 
mothur v.1.27.0

Class encapsulating Mothur galaxy tool.
Expect each invocation to include:
Here is an example call to this script with an explanation before each param :
 mothur_wrapper.py
 # name of a mothur command, this is required
 --cmd='summary.shared'
 # Galaxy output dataset list, these are output files that can be determined before the command is run
 # The items in the list are separated by commas
 # Each item contains a regex to match the output filename and a galaxy dataset filepath in which to copy the data (separated by :)
 --result='^mothur.\S+\.logfile$:'/home/galaxy/data/database/files/002/dataset_2613.dat,'^\S+\.summary$:'/home/galaxy/data/database/files/002/dataset_2614.dat
 # Galaxy output dataset extra_files_path direcotry in which to put all output files
 --outputdir='/home/galaxy/data/database/files/002/dataset_2613_files'
 # The id of one of the galaxy outputs (e.g. the mothur logfile) used for dynamic dataset generation
 #  http://wiki.galaxyproject.org/Admin/Tools/Multiple%20Output%20Files
 --datasetid='2578'
 # The galaxy directory in which to copy all output files for dynamic dataset generation
 --new_file_path='/home/galaxy/data/database/tmp'
 # specifies files to copy to the new_file_path
 # The list is separated by commas
 # Each item  conatins:   a regex pattern for matching filenames and  a galaxy datatype (separated by :)
 # The regex match.groups()[0] is used as the id name of the dataset, and must result in  unique name for each output
 --new_datasets='^\S+?\.((\S+)\.(unique|[0-9.]*)\.dist)$:lower.dist'
 ## Before version 1.18.0 Many mothur commands first required data to be read into memory using: read.otu, read.dist, or read.tree
 # This prequisite command and its params are prefixed with 'READ_'
 --READ_cmd='read.otu'
 --READ_list=/home/galaxy/data/database/files/001/dataset_1557.dat
 --READ_group='/home/galaxy/data/database/files/001/dataset_1545.dat'
 --READ_label='unique,0.07'
"""
# import pkg_resources;
import logging, os, string, sys, tempfile, glob, shutil, types, urllib, optparse, re
import shlex, subprocess
from stat import *

log = logging.getLogger( __name__ )

assert sys.version_info[:2] >= ( 2, 4 )

#debug = False
debug = True
max_processors = os.getenv('MOTHUR_MAX_PROCESSORS') if os.getenv('MOTHUR_MAX_PROCESSORS') else 8

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def __main__():
    global debug
    # tranform the logfile into html 
    # add extra file ouput
    # add object tags for svg files
    def logfile_to_html(logfile_path,htmlfile_path,tmp_input_dir_name,tmp_output_dir_name,title="Mothur Logfile"):
        if debug:  print >> sys.stdout, 'logfile_to_html %s -> %s' % (logfile_path, htmlfile_path)
        if debug:  print >> sys.stdout, 'logfile_to_html input_dir:  %s' % tmp_input_dir_name
        if debug:  print >> sys.stdout, 'logfile_to_html output_dir: %s' % tmp_output_dir_name
        txt = open(logfile_path,'r')
        in_pat = re.sub('[^a-zA-Z0-9_/]','.',tmp_input_dir_name) + '/(.+)'
        out_pat = re.sub('[^a-zA-Z0-9_/]','.',tmp_output_dir_name) + '/(.+)'
        html = open(htmlfile_path,'w')
        html.write("<html><head><title>%s</title></head>\n<body>\n<pre>\n" % title)
        try:
            for line in txt:
                if line.find('set.dir') >= 0:
                    continue
                elif line.find('put directory to ') >= 0:
                    continue
                elif line.startswith('Mothur\'s directories:') :
                    continue
                elif line.startswith('outputDir=') :
                    continue
                elif line.startswith('Type ') :
                    continue
                elif line.find(tmp_output_dir_name) >= 0:
                    # if debug:  print >> sys.stdout, 'logfile_to_html #%s#' % line
                    if line.strip().endswith('.svg'):
                        line = re.sub(out_pat,' <object id="object" type="image/svg+xml" data="\\1">\\1</object> <br><a href="\\1">\\1</a> <hr/>',line)
                    else:
                        line = re.sub(out_pat,'<a href="\\1">\\1</a>',line)
                elif line.find(tmp_input_dir_name) >= 0:
                    line = re.sub(in_pat,'\\1',line)
                html.write(line)
        except Exception, e:
            print(str(e))
            pass 
        html.write('</pre>\n</body>\n</html>\n')
        html.close()
    #Overide the processors option
    def processors_callback(option, opt, value, parser):
        val = value if value <= max_processors else max_processors
        setattr(parser.values, option.dest, val)
    #Strip out "(xx)" confidence values from a taxonomy list
    def remove_confidence_callback(option, opt, value, parser):
        # val = "'" + re.sub('\(\d+\)','',value) + "'"
        setattr(parser.values, option.dest, re.sub('\(\d+\)','',value))
    #Mothur separates items in a list of values with hyphens
    def multi_val_callback(option, opt, value, parser):
        setattr(parser.values, option.dest, '-'.join(value.split(',')))
    #Handle list arguments that may contain dataset files
    def prepare_input(value,input_dir):
        """
        Mothur requires items in a list of input values to be separated by hyphens
        Change list separator from comma to hyphen 
        If those items are files, copy them to the input dir.
        """
        if input_dir != None and value != None and isinstance(value, str): 
            vals = []
            for val in value.split(','):
                if val[0] == '/' and os.path.isfile(val):
                    (dir,fname) = os.path.split(val)
                    try:
                        os.symlink(val,os.path.join(input_dir,fname))
                    except:
                        shutil.copy(val,os.path.join(input_dir,fname))
                    vals.append(fname)
                    if debug: print >> sys.stdout, "scp %s %s" % (val, os.path.join(input_dir,fname))
                else:
                    vals.append(convert_value(val))
            return '-'.join(vals)
        return convert_value(value)
    # Ensure parameter values are in a format that mothur can handle
    def convert_value(value):
        """
        Convert parameter values to a format suitable for input to mothur
        (specifically floating point numbers supplied in scientific
        notation)
        """
        if value is None:
            # Return None
            x = None
        else:
            x = str(value)
            # Integer
            try:
                x = int(x)
            except ValueError:
                # Float
                try:
                    x = float(x)
                    if str(x).count('e'):
                        # Ugly hacks to convert scientific notation (which
                        # mothur can't handle) into decimal format
                        places = int(str(x).split('e')[1].lstrip('-'))
                        if x < 1.0:
                            x = '%.*f' % (int(places),x)
                        else:
                            x = '%*.f' % (int(places),x)
                except ValueError:
                    # Neither integer nor float
                    pass
            x = str(x)
        # Return whatever we finished up with
        return x
    #Parse Command Line and get a list of params
    # Prefix is used to differentiate between prerequisite commands: read.otu, read.dist, read.tree
    def get_params(cmd, options, input_dir, prefix=''):
        """
        Gather parameter values for the specified mothur command 'cmd',
        using the definition from the 'cmd_dict' dictionary.
        """
        if debug: print >> sys.stdout, options
        params = []  
        for opt in cmd_dict[cmd]['required']:
            if debug: print >> sys.stdout, opt
            if isinstance(opt,list): # One of these must be present
                missing = True
                for sel in opt:
                    # if debug: print >> sys.stderr, sel
                    try:
                        value = prepare_input(getattr(options,prefix+sel), input_dir)
                        if value != None:
                            params.append('%s=%s' % (sel,value))
                            missing = False
                    except Exception, e:
                        stop_err('Illegal option for cmd %s : %s %s' % (cmd,sel,e))
                if missing: 
                    stop_err('Missing a required parameter for %s, need one of %s' % (cmd,opt))
            else: # This option is required
                try:
                    value = prepare_input(getattr(options,prefix+opt), input_dir)
                    # if debug: print >> sys.stderr, value
                    if value != None:
                        params.append('%s=%s' % (opt,value))
                    else:
                       stop_err('Missing a required parameter for %s : %s' % (cmd,opt))
                except Exception, e:
                    stop_err('Illegal option %s %s' % (opt,e))
        if 'optional' in cmd_dict[cmd]:
            for opt in cmd_dict[cmd]['optional']:
                # if debug: print >> sys.stderr, opt
                try:
                    value = prepare_input(getattr(options,prefix+opt), input_dir)
                    # if debug: print >> sys.stderr, value
                    if value != None:
                        params.append('%s=%s' % (opt,value))
                except Exception, e:
                    # should distinguish between READ_ opts and cmd opts
                    # stop_err('Illegal option for %s : %s' % (cmd,opt))
                    pass
        return params
    """
    The command dict has a dict for each mothur command with required and options arguments
    The top level list of required arguments is interpreted that each value item is required and any list item represents a choice of arguments 
    This covers many, but not all of the argument dependency requirements.
    For example - read.dist  required a phylip or (column and name) argument.
    The complexity of inputs should be handled by the glaxy tool xml file.
    """
    cmd_dict = dict()
    cmd_dict['align.check'] = dict({'required' : ['fasta','map'], 'optional' : ['name','count']})
    cmd_dict['align.seqs'] = dict({'required' : ['fasta','reference',], 'optional' : ['search','ksize','align','match','mismatch','gapopen','gapextend','flip','threshold','save','processors']})
    cmd_dict['amova'] = dict({'required' : ['phylip','design'] ,  'optional' : ['alpha','iters','sets']})
    cmd_dict['anosim'] = dict({'required' : ['phylip','design'] ,  'optional' : ['alpha','iters']})
    cmd_dict['bin.seqs'] = dict({'required' : ['list','fasta'], 'optional' : ['name','label','group','count']})
    cmd_dict['bootstrap.shared'] = dict({'required' : ['shared'], 'optional' : ['calc','groups','iters','label']})
    #catchall
    cmd_dict['chimera.bellerophon'] = dict({'required' : ['fasta'], 'optional' : ['filter','correction','window','increment','processors']})
    cmd_dict['chimera.ccode'] = dict({'required' : ['fasta','reference'], 'optional' : ['filter','mask','window','numwanted','save','processors']})
    cmd_dict['chimera.check'] = dict({'required' : ['fasta','reference'], 'optional' : ['ksize','svg','name','increment','save','processors']})
    cmd_dict['chimera.perseus'] = dict({'required' : ['fasta',['name','count']], 'optional' : ['group','alpha','beta','cutoff','dereplicate']})
    cmd_dict['chimera.pintail'] = dict({'required' : ['fasta','reference'], 'optional' : ['conservation','quantile','filter','mask','window','increment','save','processors']})
    cmd_dict['chimera.slayer'] = dict({'required' : ['fasta','reference'], 'optional' : ['name','group','count','search','window','increment','match','mismatch','numwanted','parents','minsim','mincov','iters','minbs','minsnp','divergence','realign','split','blastlocation','save','processors','dereplicate']})
    cmd_dict['chimera.uchime'] = dict({'required' : ['fasta'], 'optional' : ['name','group','count','dereplicate','reference','abskew','chimealns','minh','mindiv','xn','dn','xa','chunks','minchunk','idsmoothwindow','minsmoothid','maxp','skipgaps','skipgaps2','minlen','maxlen','ucl','queryfract','processors']})
    cmd_dict['chop.seqs'] = dict({'required' : ['fasta','numbases'],  'optional' : ['countgaps','keep','short','name','group','count']})
    cmd_dict['classify.otu'] = dict({'required' : ['list','taxonomy'],'optional' : ['name','cutoff','label','group','probs','basis','reftaxonomy','count','persample']})
    cmd_dict['classify.seqs'] = dict({'required' : ['fasta','reference','taxonomy'],'optional' : ['name','search','ksize','method','match','mismatch','gapopen','gapextend','numwanted','probs','save','processors','count','relabund']})
    cmd_dict['classify.tree'] = dict({'required' : ['taxonomy','tree'],'optional' : ['name','group','cutoff']})
    #clear.memory ## not needed in galaxy framework
    cmd_dict['clearcut'] = dict({'required' : [['phylip','fasta']],'optional' : ['seed','norandom','shuffle','neighbor','expblen','expdist','ntrees','matrixout','kimura','jukes','protein','DNA']})
    cmd_dict['cluster'] = dict({'required' : [['phylip','column']] ,  'optional' : ['name','count','method','cutoff','hard','precision','sim','showabund','timing']})
    cmd_dict['cluster.classic'] = dict({'required' : ['phylip'] ,  'optional' : ['name','count','method','cutoff','hard','sim','precision']})
    cmd_dict['cluster.fragments'] = dict({'required' : ['fasta'] ,  'optional' : ['name','diffs','percent','count']})
    cmd_dict['cluster.split'] = dict({'required' : [['fasta','phylip','column']] ,  'optional' : ['name','count','method','splitmethod','taxonomy','taxlevel','showabund','cutoff','hard','large','precision','classic','timing','processors','cluster']})
    cmd_dict['collect.shared'] = dict({'required' : ['shared'], 'optional' : ['calc','label','freq','groups','all']})
    cmd_dict['collect.single'] = dict({'required' : [['list', 'sabund', 'rabund', 'shared']], 'optional' : ['calc','abund','size','label','freq']})
    cmd_dict['consensus.seqs'] = dict({'required' : ['fasta'], 'optional' : ['list','name','label','cutoff','count']})

    cmd_dict['cooccurrence'] = dict({'required' : ['shared'], 'optional' : ['iters','metric','matrixmodel','groups','label']})

    cmd_dict['corr.axes'] = dict({'required' : [['shared','relabund','metadata'],'axes'], 'optional' : ['label','groups','method','numaxes']})
    cmd_dict['count.groups'] = dict({'required' : [['group','shared','count']], 'optional' : ['accnos','groups']})
    cmd_dict['count.seqs'] = dict({'required' : ['name'], 'optional' : ['group','groups','large']})

    cmd_dict['create.database'] = dict({'required' : [['list','shared'],'repfasta','repname','constaxonomy','count'], 'optional' : ['group','label']})

    cmd_dict['degap.seqs'] = dict({'required' : ['fasta']})
    cmd_dict['deunique.seqs'] = dict({'required' : ['fasta'],  'optional' : ['name','count']})
    cmd_dict['deunique.tree'] = dict({'required' : ['tree','name'],  'optional' : []})
    cmd_dict['dist.seqs'] = dict({'required' : ['fasta'],  'optional' : ['calc','countends','output','cutoff','oldfasta','column','processors']})
    cmd_dict['dist.shared'] = dict({'required' : ['shared'], 'optional' : ['calc','label','groups','output','subsample','iters','processors']})
    cmd_dict['fastq.info'] = dict({'required' : ['fastq','format'], 'optional' : ['fasta','qfile','pacbio','oligos','bdiffs','pdiffs','tdiffs','ldiffs','sdiffs']})
    cmd_dict['filter.seqs'] = dict({'required' : ['fasta'],  'optional' : ['vertical','trump','soft','hard','processors']})
    cmd_dict['get.group'] = dict({'required' : ['shared'], 'optional' : []})
    cmd_dict['get.groups'] = dict({'required' : ['group'], 'optional' : ['groups','accnos','fasta','name','list','shared','taxonomy','design']})
    cmd_dict['get.lineage'] = dict({'required' : ['taxonomy','taxon'],'optional' : ['fasta','name','group','list','alignreport','dups']})
    cmd_dict['get.otulist'] = dict({'required' : ['list'], 'optional' : ['label','sort']})
    cmd_dict['get.oturep'] = dict({'required' : ['list',['phylip','column']], 'optional' : ['fasta','name','label','group','groups','sorted','precision','cutoff','large','weighted','count','method']})
    cmd_dict['get.otus'] = dict({'required' : ['group','list','label'], 'optional' : ['groups','accnos']})
    cmd_dict['get.rabund'] = dict({'required' : [['list','sabund']],'optional' : ['sorted','label','count']})
    cmd_dict['get.relabund'] = dict({'required' : ['shared'],'optional' : ['scale','label','groups']})
    cmd_dict['get.sabund'] = dict({'required' : [['list','rabund']],'optional' : ['label','count']})
    cmd_dict['get.seqs'] = dict({'required' : ['accnos',['fasta','qfile','name','group','list','alignreport','taxonomy']], 'optional' : ['dups','count','fastq']})
    cmd_dict['get.sharedseqs'] = dict({'required' : ['list','group'], 'optional' : ['label', 'unique', 'shared', 'output', 'fasta','shared']})
    cmd_dict['hcluster'] = dict({'required' : [['column','phylip']] , 'optional' : ['name','method','cutoff','hard','precision','sorted','showabund','timing']})
    cmd_dict['heatmap.bin'] = dict({'required' : [['list', 'sabund', 'rabund', 'shared']], 'optional' : ['label','groups','scale','sorted','numotu','fontsize']})
    cmd_dict['heatmap.sim'] = dict({'required' : [['shared','phylip','column']], 'optional' : ['calc','name','label','groups','fontsize','count']})
    cmd_dict['homova'] = dict({'required' : ['phylip','design'] ,  'optional' : ['alpha','iters','sets']})
    cmd_dict['indicator'] = dict({'required' : [['tree','design'],['shared','relabund']], 'optional' : ['groups','label','processors']})
    cmd_dict['libshuff'] = dict({'required' : ['phylip','group'],'optional' : ['groups','iters','form','sim','step','cutoff']})
    cmd_dict['list.seqs'] = dict({'required' : [['fasta','name','group','list','alignreport','taxonomy']]})
    cmd_dict['list.otulabels'] = dict({'required': [['shared','relabund','list']], 'optional': ['group','label']})
    cmd_dict['make.biom'] = dict({'required' : ['shared'] ,  'optional' : ['constaxonomy','matrixtype','groups','label','metadata','picrust','reftaxonomy']})
    cmd_dict['make.contigs'] = dict({'required' : ['ffastq','rfastq',], 'optional' : ['align','match','mismatch','gapopen','gapextend','threshold','oligos','bdiffs','pdiffs','tdiffs','processors','rindex','findex']})

    cmd_dict['make.fastq'] = dict({'required' : ['fasta','qfile'] ,  'optional' : ['format']})
    cmd_dict['make.group'] = dict({'required' : ['fasta','groups'],  'optional' : []})
    cmd_dict['make.shared'] = dict({'required' : [['list','group','biom','count']],  'optional' : ['label','groups']})
    cmd_dict['mantel'] = dict({'required' : ['phylip','phylip2'] ,  'optional' : ['method','iters']})
    cmd_dict['merge.files'] = dict({'required' : ['input','output']})
    cmd_dict['merge.groups'] = dict({'required' : [['shared','group'],'design'],  'optional' : ['groups', 'label']})
    cmd_dict['metastats'] = dict({'required' : ['shared','design'],  'optional' : ['groups', 'label','iters','threshold','sets','processors']})
    cmd_dict['nmds'] = dict({'required' : ['phylip'], 'optional' : ['axes','mindim','maxdim','iters','maxiters','epsilon']})
    cmd_dict['normalize.shared'] = dict({'required' : [['shared','relabund']], 'optional' : ['label','method','norm','groups','makerelabund']})
    cmd_dict['otu.association'] = dict({'required' : [['shared','relabund']], 'optional' : ['groups', 'label','method','metadata']})
    cmd_dict['otu.hierarchy'] = dict({'required' : ['list','label'], 'optional' : ['output']})
    cmd_dict['pairwise.seqs'] = dict({'required' : ['fasta'],  'optional' : ['align','calc','countends','output','cutoff','match','mismatch','gapopen','gapextend','processors']})
    cmd_dict['parse.list'] = dict({'required' : ['list','group'], 'optional' : ['label','count']})
    cmd_dict['parsimony'] = dict({'required' : ['tree'], 'optional' : ['group','groups','name','iters','random','processors','count']})
    cmd_dict['pca'] = dict({'required' : [['shared','relabund']], 'optional' : ['label','groups','metric']})
    cmd_dict['pcoa'] = dict({'required' : ['phylip'], 'optional' : ['metric']})

    cmd_dict['pcr.seqs'] = dict({'required' : ['fasta'], 'optional' : ['oligos','name','group','taxonomy','ecoli','start','end','nomatch','keepprimer','keepdots','processors','count','pdiffs']})

    cmd_dict['phylo.diversity'] = dict({'required' : ['tree'],'optional' : ['group','name','groups','iters','freq','scale','rarefy','collect','summary','processors','count']})
    cmd_dict['phylotype'] = dict({'required' : ['taxonomy'],'optional' : ['name','cutoff','label']})
    cmd_dict['pre.cluster'] = dict({'required' : ['fasta'],  'optional' : ['name','count','diffs','group','processors','topdown']})
    cmd_dict['rarefaction.shared'] = dict({'required' : ['shared'], 'optional' : ['calc','label','iters','groups','jumble','design','sets','groupmode','subsample','subsampleiters']})
    cmd_dict['rarefaction.single'] = dict({'required' : [['list', 'sabund', 'rabund', 'shared']], 'optional' : ['calc','abund','iters','label','freq','processors']})
    cmd_dict['remove.groups'] = dict({'required' : ['group'], 'optional' : ['groups','accnos','fasta','name','list','shared','taxonomy','design','count']})
    cmd_dict['remove.lineage'] = dict({'required' : [['taxonomy','constaxonomy'],'taxon'],'optional' : ['fasta','name','group','list','alignreport','dups','shared','list','count']})
    cmd_dict['remove.otus'] = dict({'required' : ['group','list','label'], 'optional' : ['groups','accnos']})
    cmd_dict['remove.rare'] = dict({'required' : [['list','sabund','rabund','shared'],'nseqs'], 'optional' : ['group','groups','label','bygroup','count']})
    cmd_dict['remove.seqs'] = dict({'required' : ['accnos',['fasta','qfile','name','group','list','alignreport','taxonomy']], 'optional' : ['dups','count','fastq']})
    cmd_dict['reverse.seqs'] = dict({'required' : ['fasta']})
    cmd_dict['screen.seqs'] = dict({'required' : ['fasta'],  'optional' : ['start','end','maxambig','maxhomop','minlength','maxlength','criteria','optimize','name','group','alignreport','taxonomy','processors','count','summary','contigsreport']})
    cmd_dict['sens.spec'] = dict({'required' : ['list',['column','phylip']] , 'optional' : ['label','cutoff','hard','precision']})
    cmd_dict['seq.error'] = dict({'required' : ['fasta','reference'] , 'optional' : ['name','qfile','report','ignorechimeras','threshold','processors']})
    cmd_dict['sffinfo'] = dict({'required' : [['sff','sfftxt']], 'optional' : ['oligos','bdiffs','pdiffs','tdiffs','ldiffs','sdiffs','fasta','qfile','trim','sfftxt','flow','accnos']})
    cmd_dict['shhh.flows'] = dict({'required' : [['flow']], 'optional' : ['lookup','maxiter','mindelta','cutoff','sigma','order','large','processors']})
    cmd_dict['shhh.seqs'] = dict({'required' : [['fasta','name']], 'optional' : ['group','sigma','processors']})
    cmd_dict['split.abund'] = dict({'required' : ['fasta',['name','list','count']], 'optional' : ['cutoff','group','groups','label','accnos']})
    cmd_dict['split.groups'] = dict({'required' : ['fasta','group'], 'optional' : ['name','groups','count']})
    cmd_dict['sort.seqs'] = dict({'required' : [['fasta','qfile','name','group','flow','taxonomy']], 'optional' : ['accnos','large','count']})
    cmd_dict['sub.sample'] = dict({'required' : [['fasta','list','sabund','rabund','shared']], 'optional' : ['name','group','groups','label','size','persample','count','taxonomy']})
    cmd_dict['summary.qual'] = dict({'required' : ['qfile'], 'optional' : ['name','count']})
    cmd_dict['summary.seqs'] = dict({'required' : ['fasta'], 'optional' : ['name','processors','count']})
    cmd_dict['summary.shared'] = dict({'required' : ['shared'], 'optional' : ['calc','label','groups','all','distance','processors','subsample','iters']})
    cmd_dict['summary.single'] = dict({'required' : [['list','sabund','rabund','shared']], 'optional' : ['calc','abund','size','label','groupmode','subsample','iters']})
    cmd_dict['summary.tax'] = dict({'required' : ['taxonomy'], 'optional' : ['name','group','reftaxonomy','count','relabund']})
    cmd_dict['tree.shared'] = dict({'required' : [['shared','phylip','column']], 'optional' : ['name','groups','calc','cutoff','precision','label','subsample','iters','processors','count']})
    cmd_dict['trim.flows'] = dict({'required' : ['flow'],  'optional' : ['oligos','bdiffs','pdiffs','tdiffs','ldiffs','sdiffs','minflows','maxflows','fasta','signal','noise','maxhomop','order','processors']})
    cmd_dict['trim.seqs'] = dict({'required' : ['fasta'],  'optional' : ['name','group','oligos','qfile','qaverage','qthreshold','qwindowaverage','qwindowsize','rollaverage','qstepsize','qtrim','flip','maxambig','maxhomop','minlength','maxlength','bdiffs','pdiffs','tdiffs','ldiffs','sdiffs','keepforward','allfiles','keepfirst','removelast','processors','count','checkorient','logtransform']})
    cmd_dict['unifrac.unweighted'] = dict({'required' : ['tree'], 'optional' : ['name','group','groups','iters','distance','random','root','subsample','consensus','processors','count']})
    cmd_dict['unifrac.weighted'] = dict({'required' : ['tree'], 'optional' : ['name','group','groups','iters','distance','random','root','subsample','consensus','processors','count']})
    cmd_dict['unique.seqs'] = dict({'required' : ['fasta'],  'optional' : ['name','count']})
    cmd_dict['venn'] = dict({'required' : [['list','shared']], 'optional' : ['calc','label','groups','abund','nseqs','permute','fontsize','sharedotus']})

    cmd_dict['merge.sfffiles'] = dict({'required' : ['sff','output']})

    parser = optparse.OptionParser()
    # Options for managing galaxy interaction
    parser.add_option( '--debug', dest='debug', action='store_true', default=False, help='Turn on wrapper debugging to stdout'  )
    parser.add_option( '--cmd', dest='cmd', help='The mothur command' )
    parser.add_option( '--inputdir', dest='inputdir', help='The directory in which to work' )
    parser.add_option( '--outputdir', dest='outputdir', help='The directory in which to work' )
    parser.add_option( '--tmpdir', dest='tmpdir', help='The directory in which to work' )
    parser.add_option( '--tempdefault', dest='tempdefault', help='The default directory in which to search for input' )
    parser.add_option( '--result', dest='result', help='The name pattern and destination for each output file' )
    parser.add_option( '--new_file_path', dest='new_file_path', help='The Galaxy new_file_path, dir for extra output datasets' )
    parser.add_option( '--datasetid', dest='datasetid', help='The Galaxy new_file_path, dir for extra output datasets' )
    parser.add_option( '--new_datasets', dest='new_datasets', help='The name pattern and datatype ext for ouputs' )
    # Options for prerequisite Read commands
    parser.add_option( '--READ_cmd', dest='READ_cmd', help='The mothur command' )
    # required arguments for Read.dist
    parser.add_option( '--READ_phylip', dest='READ_phylip', help='' )
    parser.add_option( '--READ_column', dest='READ_column', help='' )
    # required arguments for Read.otu
    parser.add_option( '--READ_rabund', dest='READ_rabund', help='' )
    parser.add_option( '--READ_sabund', dest='READ_sabund', help='' )
    parser.add_option( '--READ_list', dest='READ_list', help='' )
    parser.add_option( '--READ_shared', dest='READ_shared', help='' )
    parser.add_option( '--READ_relabund', dest='READ_relabund', help='' )
    # required arguments for Read.tree
    parser.add_option( '--READ_tree', dest='READ_tree', help='' )
    # optional arguments for Read cmds
    parser.add_option( '--READ_name', dest='READ_name', help='' )
    parser.add_option( '--READ_cutoff', dest='READ_cutoff',  type="float", help='' )
    parser.add_option( '--READ_hard', dest='READ_hard', help='' )
    parser.add_option( '--READ_precision', dest='READ_precision', type="int", help='' )
    parser.add_option( '--READ_sim', dest='READ_sim', help='' )
    parser.add_option( '--READ_group', dest='READ_group', help='' )
    parser.add_option( '--READ_groups', dest='READ_groups', help='' )
    parser.add_option( '--READ_ordergroup', dest='READ_ordergroup', help='')
    parser.add_option( '--READ_label', dest='READ_label', help='' )
    # Parameter specified in mothur
    parser.add_option( '--numbases', dest='numbases', type="int", help='Number of base to allow' )
    parser.add_option( '--fasta', dest='fasta', help='fasta file paths' )
    parser.add_option( '--fastq', dest='fastq', help='fastq file paths' )
    parser.add_option( '--ffastq', dest='ffastq', help='forward fastq file' )
    parser.add_option( '--rfastq', dest='rfastq', help='reverse fastq file' )
    parser.add_option( '--qfile', dest='qfile', help='Sequence read quality file (454 platform)' )
    parser.add_option( '--repfasta', dest='repfasta', help='fasta file paths' )
    parser.add_option( '--qaverage', dest='qaverage', type="int", help='Remove sequences that have an average quality below the value' )
    parser.add_option( '--qthreshold', dest='qthreshold', type="int", help='If at any point a base call in a sequence has a quality score below the value provided to the option, the sequence is terminated' )
    parser.add_option( '--qwindowaverage', dest='qwindowaverage', type="int", help='Remove sequences that have a window average quality below the value' )
    parser.add_option( '--qwindowsize', dest='qwindowsize', type="int", help='WIndow size for qwindowaverage' )
    parser.add_option( '--rollaverage', dest='rollaverage', type="int", help='Remove sequences that have a average quality below the value in a rolling window' )
    parser.add_option( '--qstepsize', dest='qstepsize', type="int", help='Distance to move a rolling window for each step' )
    parser.add_option( '--qtrim', dest='qtrim', help='For sequence below qthreshold, false to scrap file, true to trimmed and in trim file' )
    parser.add_option( '--ignorechimeras', dest='ignorechimeras', help='ignorechimeras' )
    parser.add_option( '--flip', dest='flip', help='If true, reverse complement the sequences' )
    parser.add_option( '--maxambig', dest='maxambig', type="int", help='Number of ambiguous base calls to allow' )
    parser.add_option( '--maxhomop', dest='maxhomop', type="int", help='Maximun homopolymer length allowed' )
    parser.add_option( '--minlength', dest='minlength', type="int", help='Minimun sequence length' )
    parser.add_option( '--maxlength', dest='maxlength', type="int", help='Maximun sequence length' )
    parser.add_option( '--oligos', dest='oligos', help='The oligos option takes a file that can contain the sequences of the forward and reverse primers and barcodes and their sample identifier.' )
    parser.add_option( '--bdiffs', dest='bdiffs', type="int", help='Number of barcode differences to allow' )
    parser.add_option( '--pdiffs', dest='pdiffs', type="int", help='Number of primer differences to allow' )
    parser.add_option( '--ldiffs', dest='ldiffs', type="int", help='Number of linker sequence differences to allow' )
    parser.add_option( '--sdiffs', dest='sdiffs', type="int", help='Number of spacer sequence differences to allow' )
    parser.add_option( '--tdiffs', dest='tdiffs', type="int", help='Total number of barcode and primer differences to allow' )
    parser.add_option( '--diffs', dest='diffs', type="int", help='Number of mismatched bases to allow between sequences in a group' )
    parser.add_option( '--allfiles', dest='allfiles', help='T - generate fasta and group for each barcode' )
    parser.add_option( '--keepforward', dest='keepforward', help='T - keep primer' )
    parser.add_option( '--name', dest='name', help='A file containing a 2 column table: name, and comma separated list of represetatives' )
    parser.add_option( '--repname', dest='repname', help='A file containing a 2 column table: name, and comma separated list of represetatives' )
    parser.add_option( '--accnos', dest='accnos', help='A file containing a list of names' )
    parser.add_option( '--groups', dest='groups', help='pairwise group labels' )
    parser.add_option( '--group', dest='group', help='A file containing a list of names' )
    parser.add_option( '--list', dest='list', help='A file containing a list of names' )
    parser.add_option( '--alignreport', dest='alignreport', help='A align.report file ' )
    parser.add_option( '--report', dest='report', help='' )
    parser.add_option( '--taxonomy', dest='taxonomy', help='A Taxonomy file' )
    parser.add_option( '--reftaxonomy', dest='reftaxonomy', help='A Taxonomy file' )
    parser.add_option( '--constaxonomy', dest='constaxonomy', help='The Taxonomy file output by classify.otu' )
    parser.add_option( '--taxon', dest='taxon',  help='A Taxon' )
    parser.add_option( '--taxlevel', dest='taxlevel', type="int", help='A Taxonomy level' )
    # parser.add_option( '--taxon', dest='taxon', action="callback", callback=remove_confidence_callback, help='A Taxon' )
    parser.add_option( '--candidate', dest='candidate', help=' file ' )
    parser.add_option( '--template', dest='template', help=' file ' )
    parser.add_option( '--reference', dest='reference', help=' file ' )
    parser.add_option( '--dups', dest='dups', help='if True also apply to the aliases from the names files' )
    parser.add_option( '--keep', dest='keep', help='Either front or back to specify the which end of the sequence to keep' )
    parser.add_option( '--search', dest='search', help='Method for finding the template sequence: kmer, blast, suffix' )
    parser.add_option( '--ksize', dest='ksize',  type="int", help='Size of kmers (5 - 12)' )
    parser.add_option( '--align', dest='align', help='Alignment method: needleman, blastn, gotoh' )
    parser.add_option( '--match', dest='match', help='Reward for a match, default is +1.0' )
    parser.add_option( '--mismatch', dest='mismatch', help='Penalty for a mismatch, default is -1.0' )
    parser.add_option( '--gapopen', dest='gapopen', help='Penalty for a opening, default is -2.0' )
    parser.add_option( '--gapextend', dest='gapextend',  help='Penalty for extending a gap, default is -1.0' )
    parser.add_option( '--precision', dest='precision',  help='' )
    parser.add_option( '--threshold', dest='threshold',  help='Cutoff at which an alignment is deemed bad and the reverse complement may be tried, 0.0 - 1.0 default 0.50' )
    parser.add_option( '--sim', dest='sim', help='Calculate similarity rather than distance' )
    parser.add_option( '--map', dest='map', help='File containing the secondary structure map.' )
    parser.add_option( '--label', dest='label', type='string', action="callback", callback=multi_val_callback, help='Distance levels you would like a output files created for(separated by commas or dashes)' )
    parser.add_option( '--filter', dest='filter', help='If true, a 50% soft filter will be applied' )
    parser.add_option( '--correction', dest='correction', help='If true, square root of the distances is used instead of the distance value' )
    parser.add_option( '--window', dest='window', type="int", help='Window size for searching for chimeras, default is 1/4 sequence length' )
    parser.add_option( '--increment', dest='increment', type="int", help='How far you move each window while finding chimeric sequences, default is 25' )
    parser.add_option( '--mask', dest='mask',  help='A file containing one sequence you wish to use as a mask for the your sequences' )
    parser.add_option( '--numwanted', dest='numwanted', type="int", help='How many sequences you would each query sequence compared with' )
    parser.add_option( '--start', dest='start', type="int", help='Remove sequences that start after thisposition' )
    parser.add_option( '--end', dest='end', type="int", help='Remove sequences that end before this position' )
    parser.add_option( '--criteria', dest='criteria', type="int", help='Percentage of sequences to match' )
    parser.add_option( '--optimize', dest='optimize', help='List of parameters to optimize' )
    parser.add_option( '--vertical', dest='vertical',  help='Ignore any column that only contains gap characters, "-" or "."' )
    parser.add_option( '--trump', dest='trump', help='Remove a column if the trump character is found at that position in any sequence of the alignment.' )
    parser.add_option( '--soft', dest='soft', type='int', help='Soft Mask - percentage required to retain column. (0-100)' )
    parser.add_option( '--hard', dest='hard', help='Hard Column Filter - A file should only contain one line consisting of 0 and 1 chars' )
    parser.add_option( '--calc', dest='calc', help='Calc Method - Gap Penality' )
    parser.add_option( '--count', dest='count',  help='Count file' )
    parser.add_option( '--countends', dest='countends',  help='Penalize terminal gaps' )
    parser.add_option( '--cutoff', dest='cutoff', help='Distance Cutoff threshold, discard larger distances' )
    parser.add_option( '--countgaps', dest='countgaps',  help='count gaps as bases' )
    parser.add_option( '--output', dest='output', help='Format for output' )
    parser.add_option( '--method', dest='method', help='Method to use for analysis - cluster' )
    parser.add_option( '--splitmethod', dest='splitmethod', help='Method to split a distance file - cluster.split' )
    parser.add_option( '--split', dest='split', help='Chimera split parameter, whether to detect trimeras and quadmeras' )
    parser.add_option( '--abund', dest='abund', type='int', help='Threshold for rare to Abundant OTU classification' )
    parser.add_option( '--size', dest='size', type='int', help='Size - sample size' )
    parser.add_option( '--groupmode', dest='groupmode', help='Collate groups into one result table' )
    parser.add_option( '--all', dest='all', help='Calculate for all' )
    parser.add_option( '--freq', dest='freq', type="float", help='Frequency of sequences to choose, as fraction is 0.0 - 1.0 or iteration if int > 1' )
    parser.add_option( '--iters', dest='iters', type='int', help='Iterations of randomizations' )
    parser.add_option( '--maxiter', dest='maxiter', type='int', help='Iterations' )
    parser.add_option( '--maxiters', dest='maxiters', type='int', help='Iterations of randomizations' )
    parser.add_option( '--subsample', dest='subsample', help='Number of subsample, or T to default to smallest group size' )
    parser.add_option( '--jumble', dest='jumble',  help='If false, just a collector curve across the samples' )
    parser.add_option( '--conservation', dest='conservation',  help='Template frequency information' )
    parser.add_option( '--quantile', dest='quantile',  help='Template quantile information' )
    parser.add_option( '--parents', dest='parents', type='int', help='Number of Parents to investigate' )
    parser.add_option( '--minsim', dest='minsim', type='int', help='Minimum simarity (0-100 percent)' )
    parser.add_option( '--mincov', dest='mincov', type='int', help='Minimum coverage (0-100 percent)' )
    parser.add_option( '--minbs', dest='minbs', type='int', help='Minimum bootstrap support (0-100 percent)' )
    parser.add_option( '--minsnp', dest='minsnp', type='int', help='Minimum SNPs to sample(0-100 percent)' )
    parser.add_option( '--mindim', dest='mindim', type='int', help='Minimum dimensions' )
    parser.add_option( '--maxdim', dest='maxdim', type='int', help='Maximum dimensions' )
    parser.add_option( '--percent', dest='percent', type='int', help='(0-100 percent)' )
    parser.add_option( '--divergence', dest='divergence', type='float', help='Divergence cutoff for chimera determination' )
    parser.add_option( '--sff', dest='sff',  help='Sff file' )
    parser.add_option( '--svg', dest='svg',  help='SVG' )
    parser.add_option( '--sfftxt', dest='sfftxt',  help='Generate a sff.txt file' )
    parser.add_option( '--flow', dest='flow',  help='Generate a flowgram file' )
    parser.add_option( '--minflows', dest='minflows', type='int', help='the minimum number of flows that each sequence must contain' )
    parser.add_option( '--maxflows', dest='maxflows', type='int', help='the number of flows after which all other flows should be ignored.' )
    parser.add_option( '--signal', dest='signal', type='float', help='threshold for intensity to be signal' )
    parser.add_option( '--noise', dest='noise', type='float', help='threshold for intensity to be noise' )
    parser.add_option( '--mindelta', dest='mindelta', type='float', help='threshold for determining how much change in the flowgram correction' )
    parser.add_option( '--sigma', dest='sigma', type='float', help='sigma option is used to set the dispersion of the data in the expectation-maximization' )
    parser.add_option( '--order', dest='order', help='flow order e.g. TACG' )
    parser.add_option( '--lookup', dest='lookup', help='lookup file that are needed to run shhh.seqs' )
    
    parser.add_option( '--trim', dest='trim', help='Whether sequences and quality scores are trimmed to the clipQualLeft and clipQualRight values' )
    parser.add_option( '--input', dest='input', help='' )
    parser.add_option( '--phylip', dest='phylip', help='' )
    parser.add_option( '--phylip2', dest='phylip2', help='' )
    parser.add_option( '--column', dest='column', help='' )
    parser.add_option( '--sort', dest='sort', help='specify sort order' )
    parser.add_option( '--sorted', dest='sorted', help='Input is presorted' )
    parser.add_option( '--showabund', dest='showabund', help='' )
    parser.add_option( '--short', dest='short', help='Keep sequences that are too short to chop' )
    parser.add_option( '--distance', dest='distance', help='' )
    parser.add_option( '--scale', dest='scale', help='' )
    parser.add_option( '--numotu', dest='numotu', help='' )
    parser.add_option( '--fontsize', dest='fontsize', help='' )
    parser.add_option( '--neqs', dest='neqs', help='' )
    parser.add_option( '--random', dest='random', help='' )
    parser.add_option( '--permute', dest='permute', help='' )
    parser.add_option( '--rarefy', dest='rarefy', help='' )
    parser.add_option( '--collect', dest='collect', help='' )
    parser.add_option( '--summary', dest='summary', help='' )
    parser.add_option( '--large', dest='large', help='' )
    parser.add_option( '--shuffle', dest='shuffle', help='' )
    parser.add_option( '--neighbor', dest='neighbor', help='' )
    parser.add_option( '--expblen', dest='expblen', help='' )
    parser.add_option( '--expdist', dest='expdist', help='' )
    parser.add_option( '--ntrees', dest='ntrees', help='' )
    parser.add_option( '--DNA', dest='DNA', help='' )
    parser.add_option( '--protein', dest='protein', help='' )
    parser.add_option( '--kimura', dest='kimura', help='' )
    parser.add_option( '--jukes', dest='jukes', help='' )
    parser.add_option( '--matrixout', dest='matrixout', help='' )
    parser.add_option( '--nseqs', dest='nseqs', help='' )
    parser.add_option( '--bygroup', dest='bygroup', help='' )
    parser.add_option( '--design', dest='design', help='' )
    parser.add_option( '--sets', dest='sets', help='' )
    parser.add_option( '--metric', dest='metric', help='' )
    parser.add_option( '--matrixmodel', dest='matrixmodel', help='' )
    parser.add_option( '--epsilon', dest='epsilon', help='' )
    parser.add_option( '--alpha', dest='alpha', help='' )
    parser.add_option( '--beta', dest='beta', help='' )
    parser.add_option( '--root', dest='root', help='' )
    parser.add_option( '--axes', dest='axes', help='table of name column followed by columns of axis values' )
    parser.add_option( '--numaxes', dest='numaxes', help='the number of axes' )
    parser.add_option( '--metadata', dest='metadata', help='data table with columns of floating-point values' )
    parser.add_option( '--basis', dest='basis', help='what the summary file represents' )
    parser.add_option( '--keepfirst', dest='keepfirst', help='trimming' )
    parser.add_option( '--removelast', dest='removelast', help='trimming' )
    parser.add_option( '--persample', dest='persample', help='sub sample option' )
    parser.add_option( '--timing', dest='timing', help='timing option' )
    parser.add_option( '--processors', dest='processors', type='int', action='callback', callback=processors_callback, help='Number of processors to use' )
    parser.add_option( '--dereplicate', dest='dereplicate', help='Boolean - remove chimeric sequences from all groups, default=f')
    parser.add_option( '--abskew', dest='abskew', type="float", help='Minimum abundance skew.')
    parser.add_option( '--chimealns', dest='chimealns', help='Boolean - produce a file containing multiple alignments of query sequences')
    parser.add_option( '--minh', dest='minh', type="float", help='mininum score to report chimera.')
    parser.add_option( '--mindiv', dest='mindiv', type="float", help='mininum score to report chimera.')
    parser.add_option( '--xn', dest='xn', type="float", help='weight of a no vote')
    parser.add_option( '--dn', dest='dn', type="float", help='pseudo-count prior on number of no votes')
    parser.add_option( '--xa', dest='xa', type="float", help='weight of an abstain vote')
    parser.add_option( '--chunks', dest='chunks', type="int", help='the number of chunks to extract from the query sequence')
    parser.add_option( '--minchunk', dest='minchunk', type="int", help='the minimum length of a chunk')
    parser.add_option( '--idsmoothwindow', dest='idsmoothwindow', type="int", help='ength of id smoothing window')
    parser.add_option( '--minsmoothid', dest='minsmoothid', type="float", help='minimum factional identity over smoothed window')
    parser.add_option( '--maxp', dest='maxp', type="int", help='maximum number of candidate parents to consider')
    parser.add_option( '--skipgaps', dest='skipgaps',  help='Boolean')
    parser.add_option( '--skipgaps2', dest='skipgaps2',  help='Boolean')
    parser.add_option( '--ucl', dest='ucl', help='')
    parser.add_option( '--queryfract', dest='queryfract', type="float", help='')
    parser.add_option( '--minlen', dest='minlen', type="int", help='Minimun sequence length' )
    parser.add_option( '--maxlen', dest='maxlen', type="int", help='Maximun sequence length' )
    parser.add_option( '--ecoli', dest='ecoli',  help='ecoli referance fasta' )
    parser.add_option( '--nomatch', dest='nomatch',  help='What to with non matching items' )
    parser.add_option( '--keepprimer', dest='keepprimer',  help='Whether to retain the primer' )
    parser.add_option( '--keepdots', dest='keepdots',  help='Whether to retain dots in the sequence' )
    parser.add_option( '--matrixtype', dest='matrixtype',  help='' )
    parser.add_option( '--consensus', dest='consensus',  help='boolean' )
    parser.add_option( '--biom', dest='biom',  help='biom file' )
    parser.add_option( '--classic', dest='classic',  help='boolean' )
    # include read.otu options
    parser.add_option( '--rabund', dest='rabund', help='' )
    parser.add_option( '--sabund', dest='sabund', help='' )
    parser.add_option( '--shared', dest='shared', help='' )
    parser.add_option( '--relabund', dest='relabund', help='' )
    parser.add_option( '--makerelabund', dest='makerelabund', help='Whether to convert to a relabund file' )
    parser.add_option( '--ordergroup', dest='ordergroup', help='')
    # include read.tree options
    parser.add_option( '--tree', dest='tree', help='' )

    # trim.seq options
    parser.add_option( '--checkorient', dest='checkorient', help='' )
    parser.add_option( '--logtransform', dest='logtransform', help='' )
    #fastq.info options
    parser.add_option( '--format', dest='format', help='' )
    parser.add_option( '--pacbio', dest='pacbio', help='' )
    parser.add_option( '--picrust', dest='picrust', help='' )
    parser.add_option( '--subsampleiters', dest='subsampleiters', help='' )

    (options, args) = parser.parse_args()
    """
    """
    # print >> sys.stderr, options # so will appear as blurb for file
    if debug == None and options.debug != None:
       debug = options.debug
    params = []  
    inputdir = None
    outputdir = None
    tmp_dir = None
    tempdefault = ''
    # Report exception if no command is given
    if options.cmd == None:
       stop_err('No mothur command given')
    # Read directory options
    if options.tmpdir != None:
        if not os.path.isdir(options.tmpdir):
            os.makedirs(options.tmpdir)
        tmp_dir = options.tmpdir
    else:
        if options.outputdir != None:
            if not os.path.isdir(options.outputdir):
                os.makedirs(options.outputdir)
            tmp_dir = os.path.join(options.outputdir,'tmp') 
            if not os.path.isdir(tmp_dir):
                os.makedirs(tmp_dir)
        else:
            tmp_dir = tempfile.mkdtemp()
    if options.inputdir != None:
        if not os.path.isdir(options.inputdir):
            os.makedirs(options.inputdir)
        inputdir = options.inputdir
    if options.outputdir != None:
        if not os.path.isdir(options.outputdir):
            os.makedirs(options.outputdir)
        try:
            outputdir = os.path.join(tmp_dir,'out') 
            os.symlink(options.outputdir,outputdir)
        except Exception, e:
            print >> sys.stderr, e
            outputdir = options.outputdir
    # Set directories
    if inputdir == None or not os.path.isdir(inputdir):
        inputdir = tmp_dir
    if outputdir == None or not os.path.isdir(outputdir):
        outputdir = tempfile.mkdtemp()
    if options.tempdefault != None and os.path.isdir(options.tempdefault):
        tempdefault = ', tempdefault=%s' % options.tempdefault
    # params.append('set.dir(input=%s,output=%s%s)' % (inputdir,outputdir,tempdefault))
    params.append('set.dir(output=%s%s)' % (outputdir,tempdefault))
    # Check for prerequisite read command 
    if options.READ_cmd != None:
        read_cmd_opts = ','.join(get_params(options.READ_cmd,options,inputdir,prefix='READ_'))
        params.append('%s(%s)' % (options.READ_cmd,read_cmd_opts))
    # Check for command options 
    cmd_opts = ','.join(get_params(options.cmd,options,inputdir))
    # print >> sys.stderr, cmd_opts
    # print >> sys.stderr, params # so will appear as blurb for file
    params.append('%s(%s)' % (options.cmd,cmd_opts))
    if debug: params.append('get.current()')
    try:
        # Generate the mothur commandline 
        # http://www.mothur.org/wiki/Command_line_mode
        cmdline = 'mothur "#'  + '; '.join(params) + '"'
        if debug: print >> sys.stdout, '%s' % cmdline
        if tmp_dir == None or not os.path.isdir(tmp_dir):
            tmp_dir = tempfile.mkdtemp()
        tmp_stderr_name = tempfile.NamedTemporaryFile( dir=tmp_dir,suffix='.err' ).name
        tmp_stderr = open( tmp_stderr_name, 'wb' )
        tmp_stdout_name = tempfile.NamedTemporaryFile( dir=tmp_dir,suffix='.out' ).name
        tmp_stdout = open( tmp_stdout_name, 'wb' )
        proc = subprocess.Popen( args=cmdline, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno(), stdout=tmp_stdout.fileno() )
        # proc = subprocess.Popen( args=cmdline, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE )
        returncode = proc.wait()
        if debug: print >> sys.stdout, 'returncode %d' % returncode
        tmp_stderr.close()
        # get stderr, allowing for case where it's very large
        tmp_stderr = open( tmp_stderr_name, 'rb' )
        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += tmp_stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stderr.close()
        tmp_stdout.close()
        if debug: print >> sys.stdout, 'parse %s' % tmp_stdout_name
        
        if returncode != 0:
            if returncode == -11 and options.cmd == 'seq.error':
                if debug: print >> sys.stdout, 'seq.error produced a segmentation fault but we are ignoring it.'
            else: 
                try:
                    # try to copy stdout to the logfile
                    for output in options.result.split(','):
                        # Each item has a regex pattern and a file path to a galaxy dataset
                        (pattern,path) = output.split(':')
                        if debug: print >> sys.stdout, '%s -> %s' % (pattern,path)
                        if pattern.find('\.logfile') > 0: 
                            if path != None and os.path.exists(path):
                                logfile_to_html(tmp_stdout_name,path,inputdir,outputdir,title="Mothur %s Error Logfile" % options.cmd)
                            break
                except:
                    pass
                
                raise Exception, stderr + "  Return code: " + str(returncode)
            
        stdout = ''
        # Parse stdout to provide info
        tmp_stdout = open( tmp_stdout_name, 'rb' )
        # try to find a "little" something interesting to print as info for the galaxy interface
        info = ''
        if options.cmd.startswith('chimera') and not options.cmd.endswith('check'):
            pattern = '^.*$'
            if options.cmd == 'chimera.slayer':
                # gi|11093931|MNE12|AF293003	yes
                pattern = '\S.*\tyes$' 
            elif options.cmd == 'chimera.bellerophon':
                # gi|11093939|MNB2|AF293011 is a suspected chimera at breakpoint 2195
                pattern = '\S.* suspected chimera .*'
            elif options.cmd == 'chimera.ccode':
                # gi|11093940|MNF8|AF293012 was found have at least one chimeric window.
                pattern = '\S.* chimeric window.*'
            elif options.cmd == 'chimera.pintail':
                pattern = '\S.*chimera.*Yes'
                # move new generated template .freq and .quan files to outputdir  
            chimera_count = 0
            for line in tmp_stdout:
                if re.match(pattern,line):
                    chimera_count += 1
            info += "Chimeras: %d" % chimera_count
        elif options.cmd == 'count.groups':
            fh = open(os.path.join(outputdir,'tmp.groups.count'),'w')
            for line in tmp_stdout:
                m = re.match('(.+) contains (\d+)\.',line)
                if m and len(m.groups()) == 2:
                    info += line  
                    print >> fh, "%s\t%s\n" % (m.group(1),m.group(2))
            fh.close()
        else:
            found_begin = False
            info_chars = 0
            for line in tmp_stdout:
                if re.match('mothur > ' + options.cmd + '\(.*\)', line):
                    found_begin = True
                    continue
                if line.find(outputdir) >= 0:
                    continue
                if line.startswith('**************'):
                    continue
                if re.match('^Processing.*',line):
                    continue
                if re.match('^Reading .*',line):
                    continue
                if re.match('^Merging .*',line):
                    continue
                if re.match('^DONE.*',line):
                    continue
                if re.match('.*\.\.\.\s*$',line):
                    continue
                if re.match('^\d*\s*$',line) and not line.find(' contains '):
                    continue
                # if re.match('^(unique|[0-9.]*)(\t\d+)+',line):  # abundance from cluster commands 
                if (not (options.cmd.startswith('unifrac') or options.cmd.startswith('count.groups')) 
                    and re.match('^(unique|[0-9.]+)(\t\d+)*',line)):  # abundance from cluster commands, allow unique line into info
                    continue
                if re.match('Output .*',line):
                    break
                if re.match('mothur > quit()',line):
                    break
                if found_begin and info_chars < 200:
                    info += "%s" % line
                    info_chars += len(line)
        tmp_stdout.close()
        print >> sys.stdout, info
        # Collect output files
        flist = os.listdir(outputdir)
        if debug: print >> sys.stdout, '%s' % flist
        # chimera.check can generate svg files, but they are not listed in the mothur.*.logfile, so we'll added them in here
        if options.cmd == 'chimera.check':
            svgs = []
            mothurlog = None
            for fname in flist:
                if fname.endswith('.svg'):
                    svgs.append(fname)
                elif fname.endswith('.logfile'):
                    mothurlog = fname
        # process option result first  
        # These are the known galaxy datasets listed in the --result= param
        if len(flist) > 0 and options.result:
            # items in list are separated by commas
            for output in options.result.split(','):
                # Each item has a regex pattern and a file path to a galaxy dataset
                (pattern,path) = output.split(':')
                if debug: print >> sys.stdout, '%s -> %s' % (pattern,path)
                if path == None or path == 'None':
                    continue
                found = False
                for fname in flist:
                    if debug: print >> sys.stdout, 'outdir %s match: %s' % (fname,re.match(pattern,fname))
                    if re.match(pattern,fname):
                        found = True
                        flist.remove(fname)
                        fpath = os.path.join(outputdir,fname)
                        if fname.endswith('.logfile'):  
                            # Make the logfile into html
                            logfile_to_html(fpath,path,inputdir,outputdir,title="Mothur %s Logfile" % options.cmd)
                        elif outputdir == options.outputdir:
                            # Use a hard link if outputdir is the extra_files_path, allows link from mothur logfile without copying data.
                            try:
                                if debug: print >> sys.stdout, 'link %s  %s' % (fpath, path)
                                os.link(fpath, path)
                            except:
                                if debug: print >> sys.stdout, 'copy %s  %s' % (fpath, path)
                                shutil.copy2(fpath, path)
                        else:
                            if debug: print >> sys.stdout, 'copy2 %s  %s' % (fpath, path)
                            shutil.copy2(fpath, path)
                        break
                # mothur.*.logfile may be in tmp_dir
                # chimera.pintail e.g.  generates files in the working dir that we might want to save
                if not found:
                    for fname in os.listdir(tmp_dir):
                        if debug: print >> sys.stdout, 'tmpdir %s match: %s' % (fname,re.match(pattern,fname))
                        if re.match(pattern,fname):
                            fpath = os.path.join(tmp_dir,fname)
                            if fname.endswith('.logfile'):
                                # Make the logfile into html
                                logfile_to_html(fpath,path,inputdir,outputdir,title="Mothur %s Logfile" % options.cmd)
                            else:
                                shutil.copy2(fpath, path)
                                break
        # Handle the dynamically generated galaxy datasets
        # http://bitbucket.org/galaxy/galaxy-central/wiki/ToolsMultipleOutput
        # --new_datasets=   specifies files to copy to the new_file_path
        # The list items are separated by commas
        # Each item  conatins:   a regex pattern for matching filenames and  a galaxy datatype (separated by :)
        # The regex match.groups()[0] is used as the id name of the dataset, and must result in  unique name for each output
        if options.new_datasets != None and options.new_file_path != None and options.datasetid != None:
            datasets = options.new_datasets.split(',')
            for output in options.new_datasets.split(','):
                (pattern,ext) = output.split(':');
                for fname in flist:
                    m = re.match(pattern,fname)
                    if m:
                        fpath = os.path.join(outputdir,fname)
                        if len(m.groups()) > 0:
                            # remove underscores since galaxy uses that as a field separator for dynamic datasets
                            root = m.groups()[0].replace('_','')
                        else:
                            # remove  the ext from the name if it exists, galaxy will add back later
                            # remove underscores since galaxy uses that as a field separator for dynamic datasets
                            root = re.sub('\.?'+ext+'$','',fname).replace('_','').replace('.','')
                        # filename pattern required by galaxy 
                        fn = "%s_%s_%s_%s_%s" % ( 'primary', options.datasetid, root, 'visible', ext )
                        if debug:  print >> sys.stdout, '> %s' % fpath
                        if debug:  print >> sys.stdout, '< %s' % os.path.join(options.new_file_path,fn)
                        try:
                            os.link(fpath, os.path.join(options.new_file_path,fn))
                        except:
                            shutil.copy2(fpath, os.path.join(options.new_file_path,fn))
    except Exception, e:
        msg = str(e) + stderr
        """
        if len(msg) < 50 and stdout != None:  
            # include the last line of stdout
            msg += stdout.splitlines()[-1]
        """
        stop_err( 'Error running  ' + msg)
    finally:
        # Only remove temporary directories
        # Enclose in try block, so we don't report error on stale nfs handles
        try: 
            if outputdir != options.outputdir and os.path.exists(outputdir):
                if os.path.islink(outputdir):
                    if debug:  print >> sys.stdout, 'rm outputdir %s' % outputdir
                    os.remove(outputdir)
                    if debug:  print >> sys.stdout, 'rmtree outputdir %s' % outputdir
                    shutil.rmtree(os.path.dirname(outputdir))
                else:
                    if debug:  print >> sys.stdout, 'rmtree %s' % outputdir
                    shutil.rmtree(outputdir)
            if inputdir != options.inputdir and os.path.exists(inputdir):
                if debug:  print >> sys.stdout, 'rmtree %s' % inputdir
                shutil.rmtree(inputdir)
        except:
            if debug:  print >> sys.stdout, 'rmtree failed'  
            pass
        

if __name__ == "__main__": __main__()
