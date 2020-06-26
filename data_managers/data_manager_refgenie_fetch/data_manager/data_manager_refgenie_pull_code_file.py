#!/usr/bin/env python

import logging
import refgenconf
from base64 import urlsafe_b64encode

log = logging.getLogger("tools.iuc.data_managers.data_manager_refgenie_pull")

def galaxy_code_get_refgenie_assets(refgenie_config_file):
    rgc = refgenconf.RefGenConf(refgenie_config_file)
    rval = []
    for urlname, genomes in rgc.listr().items():
        urlname_64 = urlsafe_b64encode(bytes(urlname, 'utf8')).decode('utf8')
        ul = []
        for genome, assets in genomes.items():
            al = []
            for name in assets:
                al.append({ 'name':name, 'value': '%s/%s/%s' % (urlname_64, genome, name), 'options':[], 'selected': False })
            ul.append({'name':genome, 'value':genome, 'options':al, 'selected': False })
        rval.append({'name':urlname, 'value':urlname_64, 'options':ul, 'selected': False })
    return rval
