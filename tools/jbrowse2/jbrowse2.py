#!/usr/bin/env python
# change to accumulating all configuration for config.json based on the default from the clone
import argparse
import binascii
import copy
import datetime
import hashlib
import json
import logging
import os
import shutil
import struct
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from collections import defaultdict

from Bio.Data import CodonTable
logging.basicConfig(level=logging.INFO)
log = logging.getLogger('jbrowse')
TODAY = datetime.datetime.now().strftime("%Y-%m-%d")
GALAXY_INFRASTRUCTURE_URL = None
mapped_chars = {
            '>': '__gt__',
            '<': '__lt__',
            "'": '__sq__',
            '"': '__dq__',
            '[': '__ob__',
            ']': '__cb__',
            '{': '__oc__',
            '}': '__cc__',
            '@': '__at__',
            '#': '__pd__',
            "": '__cn__'
        }

class ColorScaling(object):

    COLOR_FUNCTION_TEMPLATE = """
    function(feature, variableName, glyphObject, track) {{
        var score = {score};
        {opacity}
        return 'rgba({red}, {green}, {blue}, ' + opacity + ')';
    }}
    """

    COLOR_FUNCTION_TEMPLATE_QUAL = r"""
    function(feature, variableName, glyphObject, track) {{
        var search_up = function self(sf, attr){{
            if(sf.get(attr) !== undefined){{
                return sf.get(attr);
            }}
            if(sf.parent() === undefined) {{
                return;
            }}else{{
                return self(sf.parent(), attr);
            }}
        }};

        var search_down = function self(sf, attr){{
            if(sf.get(attr) !== undefined){{
                return sf.get(attr);
            }}
            if(sf.children() === undefined) {{
                return;
            }}else{{
                var kids = sf.children();
                for(var child_idx in kids){{
                    var x = self(kids[child_idx], attr);
                    if(x !== undefined){{
                        return x;
                    }}
                }}
                return;
            }}
        }};

        var color = ({user_spec_color} || search_up(feature, 'color') || search_down(feature, 'color') || {auto_gen_color});
        var score = (search_up(feature, 'score') || search_down(feature, 'score'));
        {opacity}
        if(score === undefined){{ opacity = 1; }}
        var result = /^#?([a-f\d]{{2}})([a-f\d]{{2}})([a-f\d]{{2}})$/i.exec(color);
        var red = parseInt(result[1], 16);
        var green = parseInt(result[2], 16);
        var blue = parseInt(result[3], 16);
        if(isNaN(opacity) || opacity < 0){{ opacity = 0; }}
        return 'rgba(' + red + ',' + green + ',' + blue + ',' + opacity + ')';
    }}
    """

    OPACITY_MATH = {
        'linear': """
            var opacity = (score - ({min})) / (({max}) - ({min}));
        """,
        'logarithmic': """
            var opacity = Math.log10(score - ({min})) / Math.log10(({max}) - ({min}));
        """,
        'blast': """
            var opacity = 0;
            if(score == 0.0) {{
                opacity = 1;
            }} else {{
                opacity = (20 - Math.log10(score)) / 180;
            }}
        """
    }

    BREWER_COLOUR_IDX = 0
    BREWER_COLOUR_SCHEMES = [
        (166, 206, 227),
        (31, 120, 180),
        (178, 223, 138),
        (51, 160, 44),
        (251, 154, 153),
        (227, 26, 28),
        (253, 191, 111),
        (255, 127, 0),
        (202, 178, 214),
        (106, 61, 154),
        (255, 255, 153),
        (177, 89, 40),
        (228, 26, 28),
        (55, 126, 184),
        (77, 175, 74),
        (152, 78, 163),
        (255, 127, 0),
    ]

    BREWER_DIVERGING_PALLETES = {
        'BrBg': ("#543005", "#003c30"),
        'PiYg': ("#8e0152", "#276419"),
        'PRGn': ("#40004b", "#00441b"),
        'PuOr': ("#7f3b08", "#2d004b"),
        'RdBu': ("#67001f", "#053061"),
        'RdGy': ("#67001f", "#1a1a1a"),
        'RdYlBu': ("#a50026", "#313695"),
        'RdYlGn': ("#a50026", "#006837"),
        'Spectral': ("#9e0142", "#5e4fa2"),
    }

    def __init__(self):
        self.brewer_colour_idx = 0

    def rgb_from_hex(self, hexstr):
        # http://stackoverflow.com/questions/4296249/how-do-i-convert-a-hex-triplet-to-an-rgb-tuple-and-back
        return struct.unpack('BBB', binascii.unhexlify(hexstr))

    def min_max_gff(self, gff_file):
        min_val = None
        max_val = None
        with open(gff_file, 'r') as handle:
            for line in handle:
                try:
                    value = float(line.split('\t')[5])
                    min_val = min(value, (min_val or value))
                    max_val = max(value, (max_val or value))

                    if value < min_val:
                        min_val = value

                    if value > max_val:
                        max_val = value
                except Exception:
                    pass
        return min_val, max_val

    def hex_from_rgb(self, r, g, b):
        return '#%02x%02x%02x' % (r, g, b)

    def _get_colours(self):
        r, g, b = self.BREWER_COLOUR_SCHEMES[self.brewer_colour_idx % len(self.BREWER_COLOUR_SCHEMES)]
        self.brewer_colour_idx += 1
        return r, g, b

    def parse_menus(self, track):
        trackConfig = {'menuTemplate': [{}, {}, {}, {}]}

        if 'menu' in track['menus']:
            menu_list = [track['menus']['menu']]
            if isinstance(track['menus']['menu'], list):
                menu_list = track['menus']['menu']

            for m in menu_list:
                tpl = {
                    'action': m['action'],
                    'label': m.get('label', '{name}'),
                    'iconClass': m.get('iconClass', 'dijitIconBookmark'),
                }
                if 'url' in m:
                    tpl['url'] = m['url']
                if 'content' in m:
                    tpl['content'] = m['content']
                if 'title' in m:
                    tpl['title'] = m['title']

                trackConfig['menuTemplate'].append(tpl)

        return trackConfig

    def parse_colours(self, track, trackFormat, gff3=None):
        # Wiggle tracks have a bicolor pallete
        trackConfig = {'style': {}}
        if trackFormat == 'wiggle':

            trackConfig['style']['pos_color'] = track['wiggle']['color_pos']
            trackConfig['style']['neg_color'] = track['wiggle']['color_neg']

            if trackConfig['style']['pos_color'] == '__auto__':
                trackConfig['style']['neg_color'] = self.hex_from_rgb(*self._get_colours())
                trackConfig['style']['pos_color'] = self.hex_from_rgb(*self._get_colours())

            # Wiggle tracks can change colour at a specified place
            bc_pivot = track['wiggle']['bicolor_pivot']
            if bc_pivot not in ('mean', 'zero'):
                # The values are either one of those two strings
                # or a number
                bc_pivot = float(bc_pivot)
            trackConfig['bicolor_pivot'] = bc_pivot
        elif 'scaling' in track:
            if track['scaling']['method'] == 'ignore':
                if track['scaling']['scheme']['color'] != '__auto__':
                    trackConfig['style']['color'] = track['scaling']['scheme']['color']
                else:
                    trackConfig['style']['color'] = self.hex_from_rgb(*self._get_colours())
            else:
                # Scored method
                algo = track['scaling']['algo']
                # linear, logarithmic, blast
                scales = track['scaling']['scales']
                # type __auto__, manual (min, max)
                scheme = track['scaling']['scheme']
                # scheme -> (type (opacity), color)
                # ==================================
                # GENE CALLS OR BLAST
                # ==================================
                if trackFormat == 'blast':
                    red, green, blue = self._get_colours()
                    color_function = self.COLOR_FUNCTION_TEMPLATE.format(**{
                        'score': "feature._parent.get('score')",
                        'opacity': self.OPACITY_MATH['blast'],
                        'red': red,
                        'green': green,
                        'blue': blue,
                    })
                    trackConfig['style']['color'] = color_function.replace('\n', '')
                elif trackFormat == 'gene_calls':
                    # Default values, based on GFF3 spec
                    min_val = 0
                    max_val = 1000
                    # Get min/max and build a scoring function since JBrowse doesn't
                    if scales['type'] == 'automatic' or scales['type'] == '__auto__':
                        min_val, max_val = self.min_max_gff(gff3)
                    else:
                        min_val = scales.get('min', 0)
                        max_val = scales.get('max', 1000)

                    if scheme['color'] == '__auto__':
                        user_color = 'undefined'
                        auto_color = "'%s'" % self.hex_from_rgb(*self._get_colours())
                    elif scheme['color'].startswith('#'):
                        user_color = "'%s'" % self.hex_from_rgb(*self.rgb_from_hex(scheme['color'][1:]))
                        auto_color = 'undefined'
                    else:
                        user_color = 'undefined'
                        auto_color = "'%s'" % self.hex_from_rgb(*self._get_colours())

                    color_function = self.COLOR_FUNCTION_TEMPLATE_QUAL.format(**{
                        'opacity': self.OPACITY_MATH[algo].format(**{'max': max_val, 'min': min_val}),
                        'user_spec_color': user_color,
                        'auto_gen_color': auto_color,
                    })

                    trackConfig['style']['color'] = color_function.replace('\n', '')
        return trackConfig


def etree_to_dict(t):
    if t is None:
        return {}

    d = {t.tag: {} if t.attrib else None}
    children = list(t)
    if children:
        dd = defaultdict(list)
        for dc in map(etree_to_dict, children):
            for k, v in dc.items():
                dd[k].append(v)
        d = {t.tag: {k: v[0] if len(v) == 1 else v for k, v in dd.items()}}
    if t.attrib:
        d[t.tag].update(('@' + k, v) for k, v in t.attrib.items())
    if t.text:
        text = t.text.strip()
        if children or t.attrib:
            if text:
                d[t.tag]['#text'] = text
        else:
            d[t.tag] = text
    return d


# score comes from feature._parent.get('score') or feature.get('score')

INSTALLED_TO = os.path.dirname(os.path.realpath(__file__))


def metadata_from_node(node):
    metadata = {}
    try:
        if len(node.findall('dataset')) != 1:
            # exit early
            return metadata
    except Exception:
        return {}

    for (key, value) in node.findall('dataset')[0].attrib.items():
        metadata['dataset_%s' % key] = value

    for (key, value) in node.findall('history')[0].attrib.items():
        metadata['history_%s' % key] = value

    for (key, value) in node.findall('metadata')[0].attrib.items():
        metadata['metadata_%s' % key] = value

    for (key, value) in node.findall('tool')[0].attrib.items():
        metadata['tool_%s' % key] = value

    # Additional Mappings applied:
    metadata['dataset_edam_format'] = '<a target="_blank" href="http://edamontology.org/{0}">{1}</a>'.format(metadata['dataset_edam_format'], metadata['dataset_file_ext'])
    metadata['history_user_email'] = '<a href="mailto:{0}">{0}</a>'.format(metadata['history_user_email'])
    metadata['hist_name'] = metadata['history_display_name']
    metadata['history_display_name'] = '<a target="_blank" href="{galaxy}/history/view/{encoded_hist_id}">{hist_name}</a>'.format(
        galaxy=GALAXY_INFRASTRUCTURE_URL,
        encoded_hist_id=metadata['history_id'],
        hist_name=metadata['history_display_name']
    )
    metadata['tool_tool'] = '<a target="_blank" href="{galaxy}/datasets/{encoded_id}/show_params">{tool_id}</a>'.format(
        galaxy=GALAXY_INFRASTRUCTURE_URL,
        encoded_id=metadata['dataset_id'],
        tool_id=metadata['tool_tool_id'],
        # tool_version=metadata['tool_tool_version'],
    )
    return metadata


class JbrowseConnector(object):

    def __init__(self, jbrowse, outdir, genomes, standalone=None, gencode=1):
        self.cs = ColorScaling()
        self.jbrowse = jbrowse
        self.outdir = outdir
        self.genome_paths = genomes
        self.standalone = standalone
        self.gencode = gencode
        self.trackIdlist = []
        self.tracksToAdd = []
        self.config_json = {}
        self.config_json_file = os.path.join(outdir, 'config.json')
        if standalone == "complete":
            self.clone_jbrowse(self.jbrowse, self.outdir)
        elif standalone == "minimal":
            self.clone_jbrowse(self.jbrowse, self.outdir, minimal=True)
        else:
            os.makedirs(self.outdir)

    def update_gencode(self):
        table = CodonTable.unambiguous_dna_by_id[int(self.gencode)]
        trackList = os.path.join(self.outdir, 'data', 'trackList.json')
        with open(trackList, 'r') as handle:
            trackListData = json.load(handle)

        trackListData['tracks'][0].update({
            'codonStarts': table.start_codons,
            'codonStops': table.stop_codons,
            'codonTable': table.forward_table,
        })

        with open(trackList, 'w') as handle:
            json.dump(trackListData, handle, indent=2)

    def subprocess_check_call(self, command, output=None):
        if output:
            log.debug('cd %s && %s >  %s', self.outdir, ' '.join(command), output)
            subprocess.check_call(command, cwd=self.outdir, stdout=output)
        else:
            log.debug('cd %s && %s', self.outdir, ' '.join(command))
            subprocess.check_call(command, cwd=self.outdir)

    def subprocess_popen(self, command):
        log.debug('cd %s && %s', self.outdir, command)
        p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, err = p.communicate()
        retcode = p.returncode
        if retcode != 0:
            log.error('cd %s && %s', self.outdir, command)
            log.error(output)
            log.error(err)
            raise RuntimeError("Command failed with exit code %s" % (retcode))

    def subprocess_check_output(self, command):
        log.debug('cd %s && %s', self.outdir, ' '.join(command))
        return subprocess.check_output(command, cwd=self.outdir)

    def _jbrowse_bin(self, command):
        return os.path.realpath(os.path.join(self.jbrowse, 'bin', command))

    def symlink_or_copy(self, src, dest):
        if 'GALAXY_JBROWSE_SYMLINKS' in os.environ and bool(os.environ['GALAXY_JBROWSE_SYMLINKS']):
            cmd = ['ln', '-s', src, dest]
        else:
            cmd = ['cp', src, dest]

        return self.subprocess_check_call(cmd)

    def process_genomes(self):
        assemblies = []
        for genome_node in self.genome_paths:
            # We only expect one input genome per run. This for loop is just
            # easier to write than the alternative / catches any possible
            # issues.
            genome_name=genome_node['meta']['dataset_dname']
            faname = genome_name + '.fa'
            fa = os.path.realpath(os.path.join(self.outdir, faname))
            shutil.copy(genome_node['path'], fa)
            faind = fa + '.fai'
            cmd = ['samtools', 'faidx', fa, "--fai-idx", faind]
            self.subprocess_check_call(cmd)
            trackDict =       {
                    "name": genome_name,
                    "sequence": {
                      "type": "ReferenceSequenceTrack",
                      "trackId": "%sReferenceSequenceTrack" % genome_name,
                      "adapter": {
                        "type": "IndexedFastaAdapter",
                        "fastaLocation": {
                          "uri":  faname,
                          "locationType": "UriLocation"
                        },
                        "faiLocation": {
                          "uri":  faname + '.fai',
                          "locationType": "UriLocation"
                        }
                      }
                    }
            }
            assemblies.append(trackDict)
            #cmd = ["jbrowse", "add-assembly", fa, "-n", genome_name, "--load", "move", "--faiLocation", faind, "--target", self.outdir,]
            #self.subprocess_check_call(cmd)
        self.config_json['assemblies'] = assemblies
        self.genome_name = genome_name
        self.genome_path = fa

    def add_default_view(self):
        cmd = ['jbrowse', 'set-default-session', '-s', self.config_json_file, '-t', ','.join(self.trackIdlist), "-n", "Default",  "--target", self.outdir,] #
        log.info('###adding default tracklist=%s' % ','.join(self.trackIdlist))
        self.subprocess_check_call(cmd)

    def write_config(self):
        with open(self.config_json_file, 'w') as fp:
            json.dump(self.config_json, fp)
        log.info('###wrote config_json=%s' % str(self.config_json))

    def read_config(self):
        with open(self.config_json_file, 'r') as fp:
            self.config_json = json.load(fp)
        log.info('### read config_json=%s' % str(self.config_json))

    def _blastxml_to_gff3(self, xml, min_gap=10):
        gff3_unrebased = tempfile.NamedTemporaryFile(delete=False)
        cmd = ['python', os.path.join(INSTALLED_TO, 'blastxml_to_gapped_gff3.py'),
               '--trim', '--trim_end', '--include_seq', '--min_gap', str(min_gap), xml]
        log.debug('cd %s && %s > %s', self.outdir, ' '.join(cmd), gff3_unrebased.name)
        subprocess.check_call(cmd, cwd=self.outdir, stdout=gff3_unrebased)
        gff3_unrebased.close()
        return gff3_unrebased.name

    def add_blastxml(self, data, trackData, blastOpts, **kwargs):
        gff3 = self._blastxml_to_gff3(data, min_gap=blastOpts['min_gap'])

        if 'parent' in blastOpts and blastOpts['parent'] != 'None':
            gff3_rebased = tempfile.NamedTemporaryFile(delete=False)
            cmd = ['python', os.path.join(INSTALLED_TO, 'gff3_rebase.py')]
            if blastOpts.get('protein', 'false') == 'true':
                cmd.append('--protein2dna')
            cmd.extend([os.path.realpath(blastOpts['parent']), gff3])
            log.debug('cd %s && %s > %s', self.outdir, ' '.join(cmd), gff3_rebased.name)
            subprocess.check_call(cmd, cwd=self.outdir, stdout=gff3_rebased)
            gff3_rebased.close()

            # Replace original gff3 file
            shutil.copy(gff3_rebased.name, gff3)
            os.unlink(gff3_rebased.name)
        url = '%s.gff3' %  trackData['label']
        dest = os.path.realpath('%s/%s' % (self.outdir,url))
        self._sort_gff(gff3, dest)
        url = url+ '.gz'
        tId =  trackData['label']
        #cmd = ['jbrowse', 'add-track', url, "--indexFile", dest +'.gz.tbi', '-n', trackData['name'], '-l', "copy" , '--trackId', trackData['label'] , "--target", self.outdir, ]
        #self.subprocess_check_call(cmd)
        #need to refine this depending on what's generated by jbrowse...
        trackDict = {
          "type": "FeatureTrack",
          "trackId": tId,
          "name": trackData['name'],
          "assemblyNames": [
            self.genome_name
          ],
          "adapter": {
            "type": "Gff3TabixAdapter",
            "gffGzLocation": {
              "locationType": "UriLocation",
              "uri": url
            },
            "index": {
              "location": {
                "locationType": "UriLocation",
                "uri": url + ".tbi"
              }
            }
          },
          "displays": [
            {
              "type": "LinearBasicDisplay",
              "displayId": "%s-LinearBasicDisplay" % tId
            },
            {
              "type": "LinearArcDisplay",
              "displayId": "%s-LinearArcDisplay" % tId
            }
          ]
        }
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)
        os.unlink(gff3)

    def add_bigwig(self, data, trackData, wiggleOpts, **kwargs):
        url = '%s.bw' %  trackData['label']
        dest = os.path.realpath('%s/%s' % (self.outdir,url))
        self.symlink_or_copy(data, dest)
        #cmd = ['jbrowse', 'add-track', dest, '-n', trackData['name'], '-l', "copy" , '--trackId', trackData['label'] , "--target", self.outdir, ]
        #self.subprocess_check_call(cmd)
        tId = trackData['label']
        trackDict =    {
              "type": "QuantitativeTrack",
              "trackId": tId,
              "name": tId,
              "assemblyNames": [
                self.genome_name,
              ],
              "adapter": {
                "type": "BigWigAdapter",
                "bigWigLocation": {
                  "locationType": "UriLocation",
                  "uri":  url
                }
              },
              "displays": [
                {
                  "type": "LinearWiggleDisplay",
                  "displayId": "%s-LinearWiggleDisplay" % tId
                }
              ]
        }
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)

    def add_bam(self, data, trackData, bamOpts, bam_index=None, **kwargs):
        url = '%s.bam' %  trackData['label']
        dest = os.path.realpath('%s/%s' % (self.outdir,url))
        self.symlink_or_copy(os.path.realpath(data), dest)
        if bam_index is not None and os.path.exists(os.path.realpath(bam_index)):
            # bai most probably made by galaxy and stored in galaxy dirs, need to copy it to dest
            self.subprocess_check_call(['cp', os.path.realpath(bam_index), dest + '.bai'])
        else:
            # Can happen in exotic condition
            # e.g. if bam imported as symlink with datatype=unsorted.bam, then datatype changed to bam
            #      => no index generated by galaxy, but there might be one next to the symlink target
            #      this trick allows to skip the bam sorting made by galaxy if already done outside
            if os.path.exists(os.path.realpath(data) + '.bai'):
                self.symlink_or_copy(os.path.realpath(data) + '.bai', dest + '.bai')
            else:
                log.warn('Could not find a bam index (.bai file) for %s', data)

        #cmd = ['jbrowse', 'add-track', dest, '-n', trackData['name'], '-l', "copy" , '--trackId', trackData['label'] , "--target", self.outdir, ]
        #self.subprocess_check_call(cmd)
        trackDict = {
          "type": "AlignmentsTrack",
          "trackId": trackData['label'],
          "name": trackData['name'],
          "assemblyNames": [
           self.genome_name
          ],
          "adapter": {
            "type": "BamAdapter",
            "bamLocation": {
              "locationType": "UriLocation",
              "uri": url
            },
            "index": {
              "location": {
                "locationType": "UriLocation",
                "uri": url + ".bai"
              }
            },
            "sequenceAdapter": {
              "type": "IndexedFastaAdapter",
              "fastaLocation": {
                "locationType": "UriLocation",
                "uri": self.genome_path
              },
              "faiLocation": {
                "locationType": "UriLocation",
                "uri": self.genome_path + ".fai"
              },
              "metadataLocation": {
                "locationType": "UriLocation",
                "uri": "/path/to/fa.metadata.yaml"
              }
            }
          }
        }
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)

        if bamOpts.get('auto_snp', 'false') == 'true':
            trackData2 = copy.copy(trackData)
            trackData2.update({
                "type": "JBrowse/View/Track/SNPCoverage",
                "key": trackData['key'] + " - SNPs/Coverage",
                "label": trackData['label'] + "_autosnp",
                "chunkSizeLimit": bamOpts.get('chunkSizeLimit', '5000000')
            })
            #self._add_track(trackData2)

    def add_vcf(self, data, trackData, vcfOpts={}, **kwargs):
        url = '%s.vcf' % trackData['label']
        dest = os.path.realpath('%s/%s' % (self.outdir,url))
        # ln?
        cmd = ['ln', '-s', data, dest]
        self.subprocess_check_call(cmd)
        cmd = ['bgzip', '-c', '>', dest]
        self.subprocess_check_call(cmd)
        cmd = ['tabix', '-p', 'vcf', dest + '.gz']
        self.subprocess_check_call(cmd)
        url = url + '.gz'
        tId = trackData['label']
        #cmd = ['jbrowse', 'add-track', dest+'.gz', '-n', trackData['name'], "--indexFile", url+'.tbi',  '-l', "copy" , '--trackId', trackData['label'] , "--target", self.outdir, ]
        #self.subprocess_check_call(cmd)
        trackDict = {
              "type": "VariantTrack",
              "trackId": tId,
              "name": trackData['name'],
              "assemblyNames": [
                self.genome_name
              ],
              "adapter": {
                "type": "VcfAdapter",
                "vcfLocation": {
                  "locationType": "UriLocation",
                  "uri":  url
                }
              },
              "displays": [
                {
                  "type": "LinearVariantDisplay",
                  "displayId": "%s-LinearVariantDisplay" % tId
                },
                {
                  "type": "ChordVariantDisplay",
                  "displayId": "%s-ChordVariantDisplay" % tId
                },
                {
                  "type": "LinearPairedArcDisplay",
                  "displayId": "%s-LinearPairedArcDisplay" % tId
                }
              ]
            }
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)

    def _sort_gff(self, data, dest):
        # Only index if not already done
        if not os.path.exists(dest+'.gz'):
            cmd = "jbrowse sort-gff %s | bgzip -c > %s.gz" % (data, dest) #"gff3sort.pl --precise '%s' | grep -v \"^$\" > '%s'"
            self.subprocess_popen(cmd)
            self.subprocess_check_call(['tabix', '-f', '-p', 'gff', dest + '.gz'])

    def _sort_bed(self, data, dest):
        # Only index if not already done
        if not os.path.exists(dest):
            cmd = ['sort', '-k1,1', '-k2,2n', data]
            with open(dest, 'w') as handle:
                self.subprocess_check_call(cmd, output=handle)

            self.subprocess_check_call(['bgzip', '-f', dest])
            self.subprocess_check_call(['tabix', '-f', '-p', 'bed', dest + '.gz'])

    def add_gff(self, data, format, trackData, gffOpts, **kwargs):
        url = '%s.gff3' % trackData['label']
        dest = os.path.realpath('%s/%s' % (self.outdir,url))
        self._sort_gff(data,dest)
        url = url + '.gz'
        tId = trackData['label']
        #cmd = ['jbrowse', 'add-track', url, "--indexFile", dest +'.gz.tbi', '-n', trackData['name'], '-l', "copy" , '--trackId', tId , "--target", self.outdir, ]
        #self.subprocess_check_call(cmd)
        trackDict =     {
              "type": "FeatureTrack",
              "trackId": tId,
              "name": trackData['name'],
              "assemblyNames": [
                self.genome_name
              ],
              "adapter": {
                "type": "Gff3TabixAdapter",
                "gffGzLocation": {
                  "locationType": "UriLocation",
                  "uri": url
                },
                "index": {
                  "location": {
                    "locationType": "UriLocation",
                    "uri": url + ".tbi"
                  }
                }
              },
              "displays": [
                {
                  "type": "LinearBasicDisplay",
                  "displayId": "%s-LinearBasicDisplay" % tId
                },
                {
                  "type": "LinearArcDisplay",
                  "displayId": "%s-LinearArcDisplay" % tId
                }
              ]
        }
        self.tracksToAdd.append(trackDict)
        self.trackIdlist.append(tId)

    def add_bed(self, data, format, trackData, gffOpts, **kwargs):
        url = '%s.bed' % trackData['label']
        dest =  os.path.realpath('%s/%s' % (self.outdir, url))
        self._sort_bed(data, dest)
        #cmd = ['jbrowse', 'add-track', dest+".gz", "--indexFile", dest +'.gz.tbi', '-n', trackData['name'], '-l', "copy" , '--trackId', trackData['label'] , "--target", self.outdir, ]
        #self.subprocess_check_call(cmd)
        tId =  trackData['label']
        trackDict = {
          "type": "FeatureTrack",
          "trackId": tId,
          "name": trackData['name'],
          "assemblyNames": [
            self.genome_name
          ],
          "adapter": {
            "type": "BedAdapter",
            "bedLocation": {
              "locationType": "UriLocation",
              "uri": url
            }
          },
          "displays": [
            {
              "type": "LinearBasicDisplay",
              "displayId": "%s-LinearBasicDisplay" % tId
            },
            {
              "type": "LinearArcDisplay",
              "displayId": "%s-LinearArcDisplay" % tId
            }
          ]
        }
        self.tracksToAdd.append(trackDict)
        log.info('added a bed - appended %s to trackstoadd =%s' % (str(trackDict),self.tracksToAdd))
        self.trackIdlist.append(tId)
        log.info('added a bed - appended %s to trackIdlist=%s' % (tId,self.trackIdlist))

    def process_annotations(self, track):
        category = track['category'].replace('__pd__date__pd__', TODAY)
        outputTrackConfig = {
            'name': track['style'].get('className', 'feature'),
            'label': track['style'].get('label', 'description'),
            'category': category,
            };

        for i, (dataset_path, dataset_ext, track_human_label, extra_metadata) in enumerate(track['trackfiles']):
            # Unsanitize labels (element_identifiers are always sanitized by Galaxy)
            for key, value in mapped_chars.items():
                track_human_label = track_human_label.replace(value, key)

            log.info('Processing %s / %s', category, track_human_label)
            outputTrackConfig['key'] = track_human_label
            # We add extra data to hash for the case of REST + SPARQL.
            if 'conf' in track and 'options' in track['conf'] and 'url' in track['conf']['options']:
                rest_url = track['conf']['options']['url']
            else:
                rest_url = ''

            # I chose to use track['category'] instead of 'category' here. This
            # is intentional. This way re-running the tool on a different date
            # will not generate different hashes and make comparison of outputs
            # much simpler.
            hashData = [str(dataset_path), track_human_label, track['category'], rest_url]
            hashData = '|'.join(hashData).encode('utf-8')
            outputTrackConfig['label'] = hashlib.md5(hashData).hexdigest() + '_%s' % i
            outputTrackConfig['metadata'] = extra_metadata
            outputTrackConfig['name'] = track_human_label

            if dataset_ext in ('gff', 'gff3'):
                self.add_gff(dataset_path, dataset_ext, outputTrackConfig,
                             track['conf']['options']['gff'])
            elif dataset_ext in ('bed', ):
                self.add_bed(dataset_path, dataset_ext, outputTrackConfig,
                             track['conf']['options']['gff'])
            elif dataset_ext == 'bigwig':
                self.add_bigwig(dataset_path, outputTrackConfig,
                                track['conf']['options']['wiggle'])
            elif dataset_ext == 'bam':
                real_indexes = track['conf']['options']['pileup']['bam_indices']['bam_index']
                if not isinstance(real_indexes, list):
                    # <bam_indices>
                    #  <bam_index>/path/to/a.bam.bai</bam_index>
                    # </bam_indices>
                    #
                    # The above will result in the 'bam_index' key containing a
                    # string. If there are two or more indices, the container
                    # becomes a list. Fun!
                    real_indexes = [real_indexes]

                self.add_bam(dataset_path, outputTrackConfig,
                             track['conf']['options']['pileup'],
                             bam_index=real_indexes[i])
            elif dataset_ext == 'blastxml':
                self.add_blastxml(dataset_path, outputTrackConfig, track['conf']['options']['blast'])
            elif dataset_ext == 'vcf':
                self.add_vcf(dataset_path, outputTrackConfig)
            else:
                log.warn('Do not know how to handle %s', dataset_ext)
            log.info('added %s - \ntrackIdlist=%s, \nconfigtracks=%s' % (dataset_ext, str(jc.trackIdlist), str(jc.config_json)))
            # Return non-human label for use in other fields


    def clone_jbrowse(self, jbrowse_dir, destination, minimal=False):
        """Clone a JBrowse directory into a destination directory.
       """
        interesting = [
                'index.html', 'static', 'version.txt'
            ]
        if minimal:
            # Should be the absolute minimum required for JBrowse to function.
            for i in interesting:
                cmd = ['cp', '-r', os.path.join(jbrowse_dir, i), destination]
                self.subprocess_check_call(cmd)
        else:
            interesting.append( 'test_data')
        for i in interesting:
                cmd = ['cp', '-r', os.path.join(jbrowse_dir, i), destination]
                self.subprocess_check_call(cmd)
        cmd = ['mkdir', '-p', os.path.join(destination, 'data', 'raw')]
        self.subprocess_check_call(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument('xml', type=argparse.FileType('r'), help='Track Configuration')

    parser.add_argument('--jbrowse', help='Folder containing a jbrowse release')
    parser.add_argument('--outdir', help='Output directory', default='out')
    parser.add_argument('--standalone', choices=['complete', 'minimal', 'data'], help='Standalone mode includes a copy of JBrowse')
    parser.add_argument('--version', '-V', action='version', version="%(prog)s 0.8.0")
    args = parser.parse_args()

    tree = ET.parse(args.xml.name)
    root = tree.getroot()

    # This should be done ASAP
    GALAXY_INFRASTRUCTURE_URL = root.find('metadata/galaxyUrl').text
    # Sometimes this comes as `localhost` without a protocol
    if not GALAXY_INFRASTRUCTURE_URL.startswith('http'):
        # so we'll prepend `http://` and hope for the best. Requests *should*
        # be GET and not POST so it should redirect OK
        GALAXY_INFRASTRUCTURE_URL = 'http://' + GALAXY_INFRASTRUCTURE_URL

    jc = JbrowseConnector(
        jbrowse=args.jbrowse,
        outdir=args.outdir,
        genomes=[
            {
                'path': os.path.realpath(x.attrib['path']),
                'meta': metadata_from_node(x.find('metadata'))
            }
            for x in root.findall('metadata/genomes/genome')
        ],
        standalone=args.standalone,
        gencode=root.find('metadata/gencode').text
    )
    jc.process_genomes()

    extra_data = {
        'visibility': {
            'default_on': [],
            'default_off': [],
            'force': [],
            'always': [],
        },
        'general': {
            'defaultLocation': root.find('metadata/general/defaultLocation').text,
            'trackPadding': int(root.find('metadata/general/trackPadding').text),
            'shareLink': root.find('metadata/general/shareLink').text,
            'aboutDescription': root.find('metadata/general/aboutDescription').text,
            'show_tracklist': root.find('metadata/general/show_tracklist').text,
            'show_nav': root.find('metadata/general/show_nav').text,
            'show_overview': root.find('metadata/general/show_overview').text,
            'show_menu': root.find('metadata/general/show_menu').text,
            'hideGenomeOptions': root.find('metadata/general/hideGenomeOptions').text,
        },
        # 'plugins': [],
        # 'plugins_python': [],
    }


    for track in root.findall('tracks/track'):
        track_conf = {}
        track_conf['trackfiles'] = []

        is_multi_bigwig = False
        try:
            if track.find('options/wiggle/multibigwig') and (track.find('options/wiggle/multibigwig').text == 'True'):
                is_multi_bigwig = True
                multi_bigwig_paths = []
        except KeyError:
            pass

        trackfiles = track.findall('files/trackFile')
        if trackfiles:
            for x in track.findall('files/trackFile'):
                if is_multi_bigwig:
                    multi_bigwig_paths.append((x.attrib['label'], os.path.realpath(x.attrib['path'])))
                else:
                    if trackfiles:
                        metadata = metadata_from_node(x.find('metadata'))

                        track_conf['trackfiles'].append((
                            os.path.realpath(x.attrib['path']),
                            x.attrib['ext'],
                            x.attrib['label'],
                            metadata
                        ))
        else:
            # For tracks without files (rest, sparql)
            track_conf['trackfiles'].append((
                '',  # N/A, no path for rest or sparql
                track.attrib['format'],
                track.find('options/label').text,
                {}
            ))

        if is_multi_bigwig:
            metadata = metadata_from_node(x.find('metadata'))

            track_conf['trackfiles'].append((
                multi_bigwig_paths,  # Passing an array of paths to represent as one track
                'bigwig_multiple',
                'MultiBigWig',  # Giving an hardcoded name for now
                {}  # No metadata for multiple bigwig
            ))

        track_conf['category'] = track.attrib['cat']
        track_conf['format'] = track.attrib['format']
        try:
            # Only pertains to gff3 + blastxml. TODO?
            track_conf['style'] = {t.tag: t.text for t in track.find('options/style')}
        except TypeError:
            track_conf['style'] = {}
            pass
        track_conf['conf'] = etree_to_dict(track.find('options'))
        jc.process_annotations(track_conf)
        print('## processed', str(track_conf),'trackIdlist',jc.trackIdlist)
    print('###done processing, trackIdlist=', jc.trackIdlist, 'config=', str(jc.config_json))
    jc.config_json['tracks'] = jc.tracksToAdd
    jc.write_config()
    jc.add_default_view()
