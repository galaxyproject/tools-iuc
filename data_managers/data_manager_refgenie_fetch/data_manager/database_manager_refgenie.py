#!/usr/bin/env python
#Dan Blankenberg

import logging
import refgenconf

log = logging.getLogger( "tools.genomespace.genomespace_exporter" )#( __name__ )

def galaxy_code_get_refgenie_folders(refgenie_config_file):
    rgc = refgenconf.RefGenConf(refgenie_config_file)
    l = rgc.listr() 
    rval = [] #defining the list
    for urlname,genomes in l.items():
        for genome, assets in genomes.items():
            al = [] #defining the list
            for name in assets:
                al.append({ 'name':name, 'value': '%s/%s' % (genome,name), 'options':[], 'selected': False })
            rval.append({'name':genome, 'value': genome, 'options':al, 'selected': False })
    # cur_options.append( { 'name':directory.get( 'name' ), 'value': directory.get( 'path'), 'options':[], 'selected': selected  } )
    return rval
