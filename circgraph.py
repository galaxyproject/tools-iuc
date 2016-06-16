#!/usr/bin/env python
import itertools
import sys
import copy
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

class DataInterpreter():

    def __init__(self, xmlroot, files_dict=[]):
        self.xmlroot = xmlroot
        self.object_list = []
        self.base_path = 'output'

        self._create_object_list(xmlroot)
        self.files_dict = self._get_files_from_xml()

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

            if file_type == "karyotype":
                self._make_karyotype(file_path)
            elif file_type == "scatter":
                if ext in ('.gff', '.gff3'):
                    self.parse_gff3(file_type)
                elif ext in ('.wig', '.bigWig', '.bw'):
                    self.parse_Wig(file_type)
            elif ext == ".bed":
                result = self.parse_bed(self.files_dict[file_type])
            elif ext in (".bw", ".bigWig", '.wig'):
                result = self.parse_Wig(self.files_dict[file_type])
            elif ext in (".gff", ".gff3"):
                result = self.parse_gff3(file_type)
            elif ext == ".fa":
                result = self.parse_fasta(self.files_dict[file_type])

            if result:
                replacement_files_data[file_type] = result


    def _create_object_list(self, xmlroot):
        current_obj = None
        current_sub = None
        subobjind= 0
        recognized_objects = {
            # 'break_style': Break_Style(),
            # 'highlight': Highlight(),
            'ideogram': Ideogram(),
            'image': Image(),
            # 'karyotype': Karyotype(),
            # 'link': Link(),
            # 'pairwise': Pairwise(),
            'plot': Plot(),
            # 'rule': Rule(),
            # 'tick': Tick(),
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

    def _make_karyotype(self, fasta):
        self.karyotype = 'karyotype.txt'

        with open(os.path.join(self.base_path, self.karyotype), 'w') as output:
            data = []
            for idx, seq in enumerate(SeqIO.parse(fasta, 'fasta')):
                data.append(' '.join(map(str, [
                    'chr',
                    '-',
                    seq.id,
                    idx,
                    0,
                    len(seq),
                    COLORS_LIST[idx % len(COLORS_LIST)]
                ])))
            output.write('\n'.join(data))

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
        self.seq_dict = {}
        seqname = ''
        seqres = ''
        with open(f) as input_handle:
            for line in input_handle.readlines():
                if line[0] == '>':
                    if seqname != '':
                        self.seq_dict[seqname] = seqres
                        seqres = ''
                    seqname = line[1:].strip()
                else:
                    seqres += line.strip()
        self.seq_dict[seqname] = seqres
        return

    def parse_Wig(self, sourcekey):
        wigfile = self.files_dict[sourcekey]
        features_dict = {}
        locidict = {}
        path = os.path.dirname(wigfile)
        name, ext = os.path.splitext(wigfile)
        newfile = name + '.wig'

        if ext != '.wig':
            subprocess.check_call([
                'bigWigToWig',
                wigFile,
                newfile
            ])
        else:
            shutil.copy(wigfile, self.base_path)

        with open(newfile, 'r') as handle:
            i = 0
            currentchrom = ''
            currentspan = 0
            start = 1
            step = 1
            m = 1
            mode = ''
            locidict = {}
            tmpdict = {}
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
                    asdf = start + m * step
                    key = '%s_%s-%s' % (currentchrom, asdf, asdf + currentspan)
                    tmpdict[key] = {
                        'chromosome': currentchrom,
                        'start': asdf,
                        'end': asdf + currentspan,
                        'val': l[0]
                    }
                    m+=1
                i+=1

            for key in tmpdict:
                if tmpdict[key]['chromosome'] not in locidict.keys():
                    locidict[tmpdict[key]['chromosome']] = {}

                locidict[tmpdict[key]['chromosome']]['key'] = tmpdict[key]

        if sourcekey == 'scatter':
            with open('scatter.txt', 'w') as handle:
                for genome in locidict:
                    for feature in locidict[genome]:
                        handle.write(' '.join(map(str, [
                            locidict[genome][feature]['chromosome'],
                            locidict[genome][feature]['start'],
                            locidict[genome][feature]['end'],
                            locidict[genome][feature]['val'],
                        ])) + '\n')

        self.files_dict[sourcekey] = 'scatter.txt'

    def parse_bed(self, f):
        bed_standard_fields = ['chromosome', 'start', 'end',
                       'name', 'score', 'strand',
                       'thickstart', 'thickend', 'rgb',
                       'blockcount', 'blocksizes', 'blockstarts']
        features_dict = {}
        with open(f) as tmp:
            for l in tmp.readlines():
                data = l.strip().split()
                tmpdict = dict(zip(bed_standard_fields, data))
                identifier = tmpdict['chromosome']+'_'+tmpdict['start']+'-'+tmpdict['end']
                features_dict[identifier] = tmpdict
        self.files_dict[f] = features_dict

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

class Karyotype(TopLevelObj):
    def __init__(self):
        TopLevelObj.__init__(self)
        self.__name__ = 'karyotype'

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
    def __init__(self):
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

class CircosPlot():
    def __init__(self, basefilename, obj_list, karyotype, files_dict):
        self.base_path = 'output'
        self.files_dict = files_dict
        self.karyotype = karyotype
        self.basename = basefilename
        self.data_dict = {}
        self.last_filled_radius = 1.1
        self.master_struct = {'ideograms':{}} #Will be loop-able for jinja2
        self.plots = ['<plots>']
        self.track_dict = {} #Temporary -- will be obviated by XML at later stage.
        self.add_ideogram()
        self.files_list = []
        self.add_top_lvl_params()
        self.obj_list = obj_list
        self.boilerplate = [
            '<<include colors_fonts_patterns.conf>>',
            '<<include housekeeping.conf>>'
        ]

    def write_conf(self):
        conf_data = []

        for elem in self.boilerplate:
            conf_data.append(elem)

        current_obj = None
        prev_obj = None

        conf_data.append('karyotype = %s' % self.karyotype)
        for obj in self.obj_list:
            if obj.__name__ not in ('ideogram', 'rule', 'image', 'plot', 'zoom', 'highlight', 'tick', 'link', 'karyotype'):
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
                conf_data.append('<spacing>\ndefault = 0.25r\n</spacing>')

            for attr in dir(obj):
                if attr == 'boilerplate':
                    conf_data.append(obj.boilerplate.strip())

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

        with open(os.path.join(self.base_path, self.basename+'.conf'), 'w') as handle:
            for line in conf_data:
                handle.write(line + '\n')


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

    def append_data(self, interpreter):
        #Add interpreter data to ConfWriter dict. Depends on what kind of data structure we decide on...
        self.raw_dict = interpreter.files_dict
        self.seq_dict = interpreter.seq_dict
        self._make_karyotype()

    def add_ideogram(self, spacing_opt='0.005r', radius_opt='0.5r', thickness_opt='20p', fill_opt='yes', color_opt=[]):
        l = len(self.master_struct['ideograms']) + 1
        ideogramid = 'ideogram-'+str(l)
        self.master_struct['ideograms'][ideogramid] = {}
        self.add_spacing(spacing_opt, ideogramid)
        block = {'radius':radius_opt,
             'thickness':thickness_opt,
             'fill':fill_opt}
        self.master_struct['ideograms'][ideogramid].update(block)

    def add_spacing(self, spacing_opt, ideogramid):
        self.master_struct['ideograms'][ideogramid].update({'spacing':{'default':spacing_opt}})

    def create_base_conf_template(self):
        E = Environment(loader=PackageLoader('dummy', 'Templates'))
        T = E.get_template('template.conf')
        with open(self.basename + '.conf', 'w') as f:
            rendering = T.render(master_struct = self.master_struct, filename = self.basename + '.png', karyotype_filename = self.karyotype_filename)
            f.write(rendering)

    def _make_karyotype(self, colors_list=[]):
        if colors_list == []:
            colors_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
                       'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                       'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18'
                       'chr19', 'chr20', 'chr21', 'chr22', 'chrx', 'chry'
                    ] #Temporary -- will update to include more diverse list of colors, Brewer...
        self.karyotype_filename = self.basename + '-karyotype.txt'
        i = 0
        with open(self.karyotype_filename, 'w') as karyotype:
            for element in self.seq_dict:
                karyotype.write('chr - '+element+' '+element+' 0 '+str(len(self.seq_dict[element]))+' '+colors_list[i]+'\n')
                if i < len(colors_list)-1:
                    i+=1
                else:
                    raise ValueError('Colors list exhausted')

    def add_four_elem_track(self, type_='four-element', filename=''):
        #Generic for heatmap, scatterplot, and histogram
        if 'plots' not in self.master_struct:
            self.master_struct['plots'] = {}
        l = len(self.master_struct['plots']) + 1
        plotname = 'plot-'+type_+'-'+str(l)
        self.master_struct['plots'][plotname] = {}
        if type_ not in self.track_dict:
            key = type_
        else:
            i = 1
            while type_ +'-'+str(i) in self.track_dict:
                i+=1
            key = type_+'-'+str(i)
        if filename == '':
            f = self.basename+'-'+key+'.txt'
        else:
            f = filename
        self.track_dict[key] = ''
        self._four_element_processing(f, key)
        self._add_plot(type_, f, plotname)

    def _four_element_processing(self, file_, key):
        #Heatmaps, scatterplots, and histograms share the same 4-element format...
        with open(file_, 'w') as handle:
            for f in self.raw_dict:
                if os.path.splitext(f)[1] in ('.wig', '.bw', '.bigWig'):
                    for chromosome in self.raw_dict[f]:
                        for feature in self.raw_dict[f][chromosome]:
                            datastring = chromosome+' '+str(self.raw_dict[f][chromosome][feature]['start'])+' '+str(self.raw_dict[f][chromosome][feature]['end'])+' '+str(self.raw_dict[f][chromosome][feature]['val'])+'\n'
                            handle.write(datastring)

    def update_master_conf(self):
        self.plots.append('</plots>')
        with open(self.mainconf, 'a') as conf:
            for element in self.plots:
                conf.write(element+'\n')

    def _add_plot(self, block_type, filename, plotname, extend_bin='no',
              max_ = '1.0', min_ = '0.0', glyph='rectangle',
              glyph_size='8', color_list='spectral-9-div',
              stroke_color='dred', thickness='1', color='red', log_base_opt='1'):
        r1 = str(self.last_filled_radius + 0.100) + 'r'
        r2 = str(self.last_filled_radius + 0.200) + 'r'
        self.last_filled_radius += 0.20
        if block_type == 'hist':
            block = {'type':'histogram',
                 'file':filename,
                 'r0':r1,
                 'r1':r2,
                 'extend_bin':extend_bin}

        elif block_type == 'heat':
            block = {'show':'yes',
                 'type':'heatmap',
                 'file':filename,
                 'color':color_list,
                 'stroke_thickness':thickness,
                 'r0':r1,
                 'r1':r2,
                 'scale_log_base':log_base_opt
                 }
        elif block_type == 'scatter':
            block = {'show':'yes',
                 'type':'scatter',
                 'file':filename,
                 'r0':r1,
                 'r1':r2,
                 'max':max_,
                 'min':min_,
                 'glyph':glyph,
                 'glyph_size':glyph_size,
                 'color':color,
                 'stroke_color':stroke_color,
                 'stroke_thickness':thickness,
                 }
        else:
            raise ValueError('Unsupported plot type.')
        self.master_struct['plots'][plotname] = block

    def add_links(self):
        pass

if __name__ == "__main__":
    xml_file = sys.argv[1]
    tree = xml.etree.ElementTree.parse(xml_file)
    root = tree.getroot()

    di = DataInterpreter(root, [])
    di.files_dict = {
        'karyotype':'./test-data/miro.fa',
        'scatter':'./test-data/miro.wig'
    }
    di._parse_files()

    plotter = CircosPlot('circos', di.object_list, di.karyotype, di.files_dict)
    plotter.write_conf()
