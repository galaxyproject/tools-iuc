#!/usr/bin/env python
import itertools
import sys
import os
import re
import xml
import brewer2mpl
import shutil
import subprocess
import argparse
from jinja2 import Template, Environment, PackageLoader
from BCBio import GFF
from subprocess import call
from Bio.Seq import Seq
from Bio import SeqIO
import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()

COLORS_LIST = ['red', 'green', 'blue', 'purple', 'orange', 'yellow', 'black',
    'grey', 'white', 'pred', 'pgreen', 'pblue', 'ppurple', 'pyellow', 'porange']

BED_STANDARD_FIELDS = [
    'chromosome', 'start', 'end', 'name', 'score', 'strand', 'thickstart',
    'thickend', 'rgb', 'blockcount', 'blocksizes', 'blockstarts'
]


class CircosPlotter:

    def __init__(self, xmlroot, files_dict=None, output=None):
        self.xmlroot = xmlroot
        self.object_list = []

        self.directory = output
        self._create_object_list(xmlroot)
        if files_dict:
            self.files_dict = files_dict
        else:
            self.files_dict = self._get_files_from_xml()

        self.basename = 'circos'
        self.data_dict = {}
        self.last_filled_radius = 1.1
        self.master_struct = {'ideograms':{}} #Will be loop-able for jinja2
        self.plots = ['<plots>']
        self.track_dict = {} #Temporary -- will be obviated by XML at later stage.
        self.add_ideogram()
        self.files_list = []
        self.add_top_lvl_params()
        self.obj_list = self.object_list
        self._parse_files()

    def write_conf(self):
        conf_data = [
            '<<include colors_fonts_patterns.conf>>',
            '<<include housekeeping.conf>>'
        ]

        current_obj = None
        prev_obj = None

        conf_data.append('karyotype = karyotype.txt')
        for obj in self.obj_list:
            if obj.__name__ not in ('ideogram', 'rule', 'image', 'plot', 'zoom', 'highlight', 'tick', 'link'):
                raise Exception()

            if current_obj is not None:
                conf_data.append('</%s>' % current_obj)

            prev_obj = current_obj
            current_obj = obj.__name__

            if current_obj == 'plot':
                if prev_obj != current_obj:
                    conf_data.append('<plots>')
                else:
                    conf_data.append('</plots>')

            conf_data.append('<%s>' % current_obj)

            if current_obj == 'ideogram':
                conf_data += ['<spacing>', 'default = 0.25r', '</spacing>']

            for attr in dir(obj):
                if attr == 'boilerplate':
                    conf_data.append(obj.boilerplate.strip())
                    continue
                elif attr.startswith('__') or \
                    attr in ('llo',) or \
                    getattr(obj, attr) is None or \
                    hasattr(attr, '__call__'):
                        continue

                try:
                    # if getattr(obj, attr).endswith('__'):
                        # val = str(getattr(obj, attr)).strip()[-6:]
                        # h = [val[0:2], val[2:4], val[4:6]]
                        # rgb = []
                        # for num in h:
                            # rgb.append(str(int(str(num), 16)))
                        # s = ', '
                        # val = s.join(rgb)
                    # else:
                    val = str(getattr(obj, attr)).strip()

                    if attr[0:4] == obj.__name__[0:4]:
                        if attr[5:] in ['thickness', 'radius', 'r0', 'r1']:
                            conf_data.append('%s = %sr' % (attr[5:], val.strip()))
                        else:
                            conf_data.append('%s = %s' % (attr[5:], val.strip()))
                    else:
                        conf_data.append('%s = %s' % (attr, val.strip()))
                except Exception:
                    pass

        conf_data.append('</%s>' % current_obj)
        if current_obj == 'plot':
            conf_data.append('</plots>')

        with open(os.path.join(self.directory, self.basename + '.conf'), 'w') as handle:
            indentation_level = 0
            for line in conf_data:
                if line.startswith('<<'):
                    pass
                elif line.startswith('</'):
                    indentation_level -= 1
                elif line.startswith('<'):
                    if indentation_level == 0:
                        handle.write('\n')

                handle.write(('  ' * indentation_level) + line + '\n')

                if line.startswith('<<'):
                    pass
                elif line.startswith('</'):
                    # blank line after closing group
                    pass
                elif line.startswith('<'):
                    indentation_level += 1

    def add_top_lvl_params(self):
        self.anglestep = None
        self.beziersamples = None
        self.chromosomes_display_default = None
        self.chromosomes = None
        self.chromosomes_breaks = None
        self.chromosomes_order = None
        self.chromosomes_radius = None
        self.chromosomes_reverse = None
        self.chromosomes_units = None
        self.debug = None
        self.imagemap = None
        self.minslicestep = None
        self.show_ticks = None
        self.show_tick_labels = None
        self.units_nounit = None
        self.units_ok = None
        self.warnings = None

    def add_ideogram(self, spacing_opt='0.005r', radius_opt='0.5r',
            thickness_opt='20p', fill_opt='yes'):

        currentLength = len(self.master_struct['ideograms']) + 1
        ideogramid = 'ideogram-%s' % currentLength

        block = {
            'spacing': {
                'default': spacing_opt
            },
            'radius':radius_opt,
            'thickness':thickness_opt,
            'fill':fill_opt,
        }
        self.master_struct['ideograms'][ideogramid] = block

    def _parse_files(self):
        replacement_files_data = {}
        for (file_type, file_path) in self.files_dict.iteritems():
            log.debug("_parse_files %s %s", file_type, file_path)
            #TODO Each parse method needs to return a txt file compatible with Circos.
            name, ext = os.path.splitext(file_path)

            fileinfo = {
                'name': name,
                'ext': ext,
                'type': file_type,
                'path': file_path,
            }

            result = None

            print fileinfo
            if file_type == "scatter":
                if ext in ('.gff', '.gff3'):
                    self.parse_gff3(file_type)
                elif ext in ('.wig', '.bigWig', '.bw'):
                    result = self.parse_Wig(fileinfo)
            elif ext == ".bed":
                result = self.parse_bed(fileinfo)
            elif ext in (".bw", ".bigWig", '.wig'):
                result = self.parse_Wig(fileinfo)
            elif ext in (".gff", ".gff3"):
                result = self.parse_gff3(file_type)
            elif ext == ".fa":
                result = self.parse_fasta(self.files_dict[file_type])
                self.seq_dict = result

            if result:
                replacement_files_data[file_type] = result

        self.files_dict = replacement_files_data

    def _create_object_list(self, xmlroot):
        current_obj = None
        current_sub = None
        subobjind= 0
        recognized_objects = {
            # 'break_style': Break_Style(),
            # 'highlight': Highlight(),
            'ideogram': Ideogram(),
            'image': Image(),
            # 'link': Link(),
            # 'pairwise': Pairwise(),
            'plot': Plot(),
            # 'rule': Rule(),
            'tick': Tick(),
            # 'zoom': Zoom(),
            # 'spacing': Spacing(),
        }
        for element in xmlroot.iter():
            if element.tag in recognized_objects:
                if current_obj is not None:
                    self.object_list.append(current_obj)

                current_obj = recognized_objects[element.tag]
                subobjind = 0
            elif current_obj is not None:
                if subobjind == 0:
                    setattr(current_obj, element.tag, element.text)
                elif subobjind == 1:
                    setattr(current_sub, element.tag, element.text)
                    current_obj.llo.append(current_sub)
        self.object_list.append(current_obj)

    def _get_files_from_xml(self):
        #TODO When the time comes to create the final Galaxy tool...
        pass

    def parse_gff3(self, key):
        #FIXME Seems there's no nontrivial way to move gff3 file info into Circos... For now...
        raise Exception()
        with open(self.files_dict[key], 'r') as input_handle:
            for rec in GFF.parse(input_handle):
                flat_list = []
                features_dict = {}
                tmpdict = {
                    'id': rec.id,
                    'annotations': rec.annotations,
                    'dbxrefs': rec.dbxrefs,
                    'description': rec.description,
                    'letter_annotations': rec.letter_annotations,
                    'name': rec.name,
                    'seq': rec.seq,
                }
                for obj in rec.features:
                    if obj.sub_features == []:
                        flat_list.append(obj)
                    else:
                        flat_list = self._process_obj_subfeatures(obj, flat_list)
                tmpdict['features'] = flat_list
                features_dict[rec.id] = tmpdict
        if key == 'scatter':
            filename = str(key)+'.txt'
            g = open(filename, 'w')
            g.close()
        elif key == 'heat':
            pass
        else:
            pass
        self.files_dict[key] = g

    def _process_obj_subfeatures(self, obj, flist):
        flist.append(obj)
        l = obj.sub_features
        for element in l:
            if l != []:
                self._process_obj_subfeatures(element, flist)
            else:
                flist.append(element)
        return flist

    def parse_fasta(self, f):
        sequences = {}
        for seq in SeqIO.parse(f, 'fasta'):
            sequences[seq] = str(seq.seq)
        return sequences

    def parse_Wig(self, fileinfo):
        wigfile = fileinfo['path']
        # TODO: replace with tempfile


        if fileinfo['ext'] != '.wig':
            newfile = fileinfo['path'] + '.wig'
            subprocess.check_call([
                'bigWigToWig',
                wigFile,
                newfile
            ])
        else:
            newfile = wigfile

        i = 0
        currentchrom = ''
        currentspan = 0
        start = 1
        step = 1
        m = 1
        mode = ''
        locidict = {}
        tmpdict = {}

        with open(newfile, 'r') as handle:
            for line in handle:
                l = line.split()
                if l[0] == 'variableStep':
                    mode = 'variable'
                    currentchrom = l[1].split('=')[1]
                    try:
                        currentspan = int(l[2].split('=')[1])
                    except IndexError:
                        currentspan = 0
                elif l[0] == 'fixedStep':
                    mode = 'fixed'
                    m = 1
                    currentchrom = l[1].split('=')[1]
                    start = int(l[2].split()[1])
                    step = int(l[3].split()[1])
                    try:
                        currentspan = int(l[4].split('=')[1])
                    except IndexError:
                        currentspan = 0
                elif mode == 'variable':
                    l0 = int(l[0])
                    key = '%s_%s-%s' % (currentchrom, l[0], l0 + currentspan)
                    tmpdict[key] = {
                        'chromosome': currentchrom,
                        'start': l0,
                        'end': l0 + currentspan,
                        'val': l[1]
                    }
                elif mode == 'fixed':
                    fdsa = start + m * step
                    key = '%s_%s-%s' % (currentchrom, fdsa, fdsa + currentspan)
                    tmpdict[key] = {
                        'chromosome': currentchrom,
                        'start': fdsa,
                        'end': fdsa + currentspan,
                        'val': l[0]
                    }
                    m+=1
                i+=1

            for key in tmpdict:
                if tmpdict[key]['chromosome'] not in locidict.keys():
                    locidict[tmpdict[key]['chromosome']] = {}

                locidict[tmpdict[key]['chromosome']][key] = tmpdict[key]

        if fileinfo['type'] == 'scatter':
            with open(os.path.join(self.directory, 'scatter.txt'), 'w') as handle:
                for genome in locidict:
                    for feature in locidict[genome]:
                        handle.write(' '.join(map(str, [
                            locidict[genome][feature]['chromosome'],
                            locidict[genome][feature]['start'],
                            locidict[genome][feature]['end'],
                            locidict[genome][feature]['val'],
                        ])) + '\n')

        return 'scatter.txt'

    def parse_bed(self, fileinfo):
        features_dict = {}
        with open(fileinfo['path'], 'r') as tmp:
            for line in tmp:
                data = line.strip().split()
                tmpdict = dict(zip(BED_STANDARD_FIELDS, data))
                identifier = '{chromosome}_{start}-{end}'.format(**tmpdict)
                features_dict[identifier] = tmpdict
        return features_dict



class CircosObj():
    def __init__(self):
        pass


class TopLevelObj(CircosObj):
    def __init__(self):
        self.llo = []
        pass


class Highlight(TopLevelObj):
    def __init__(self):
        TopLevelObj.__init__(self)
        self.__name__ = 'highlight'


class Ideogram(TopLevelObj):
    def __init__(self):
        TopLevelObj.__init__(self)
        self.__name__ = 'ideogram'


class Plot(TopLevelObj):
    def __init__(self):
        TopLevelObj.__init__(self)
        self.__name__ = 'plot'


class Zoom(TopLevelObj):
    def __init__(self):
        TopLevelObj.__init__(self)
        self.__name__ = 'zoom'


class Link(TopLevelObj):
    def __init__(self):
        TopLevelObj.__init__(self)
        self.__name__ = 'link'


class Tick(TopLevelObj):
    def __init__(self):
        TopLevelObj.__init__(self)
        self.__name__ = 'tick'


class Image(TopLevelObj):
    def __init__(self):
        TopLevelObj.__init__(self)
        self.__name__ = 'image'
        self.boilerplate = '<<include etc/image.conf>>'


class LowLevelObj(CircosObj):
    pass


class Pairwise(LowLevelObj):
    def __init__(self):
        self.chromosome_pair = ['', '']
        self.__name__ = 'pairwise'


class Spacing(LowLevelObj):
    def __init__(self):
        self.__name__ = 'spacing'


class Break_Style(LowLevelObj):
    def __init__(self):
        self.__name__ = 'break_style'


class Rule(LowLevelObj):
    def __init__(self):
        self.__name__ = 'rule'

    def associate_parent(self, parent):
        self.Parent = parent #Not to be added to conf file.




if __name__ == "__main__":
    xml_file = sys.argv[1]
    tree = xml.etree.ElementTree.parse(xml_file)
    root = tree.getroot()

    di = CircosPlotter(root, {
        'scatter':'./test-data/miro.wig'
    }, output='output')
    di.write_conf()
