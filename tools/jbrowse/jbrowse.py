#!/usr/bin/env python
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
import tempfile
import xml.etree.ElementTree as ET
from collections import defaultdict

from Bio.Data import CodonTable
logging.basicConfig(level=logging.INFO)
log = logging.getLogger('jbrowse')
TODAY = datetime.datetime.now().strftime("%Y-%m-%d")
GALAXY_INFRASTRUCTURE_URL = None


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
        self.tracksToIndex = []

        if standalone == "complete":
            self.clone_jbrowse(self.jbrowse, self.outdir)
        elif standalone == "minimal":
            self.clone_jbrowse(self.jbrowse, self.outdir, minimal=True)
        else:
            try:
                os.makedirs(self.outdir)
            except OSError:
                # Ignore if the folder exists
                pass

            try:
                os.makedirs(os.path.join(self.outdir, 'data', 'raw'))
            except OSError:
                # Ignore if the folder exists
                pass

        self.process_genomes()
        self.update_gencode()

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
        for genome_node in self.genome_paths:
            # We only expect one input genome per run. This for loop is just
            # easier to write than the alternative / catches any possible
            # issues.

            # Copy the file in workdir, prepare-refseqs.pl will copy it to jbrowse's data dir
            local_genome = os.path.realpath('./genome.fasta')
            shutil.copy(genome_node['path'], local_genome)

            cmd = ['samtools', 'faidx', local_genome]
            self.subprocess_check_call(cmd)

            self.subprocess_check_call([
                'perl', self._jbrowse_bin('prepare-refseqs.pl'),
                '--trackConfig', json.dumps({'metadata': genome_node['meta']}),
                '--indexed_fasta', os.path.realpath(local_genome)])

            os.unlink(local_genome)
            os.unlink(local_genome + '.fai')

    def generate_names(self):
        # Generate names
        args = [
            'perl', self._jbrowse_bin('generate-names.pl'),
            '--hashBits', '16'
        ]

        tracks = ','.join(self.tracksToIndex)
        if tracks:
            args += ['--tracks', tracks]
        else:
            # No tracks to index, index only the refseq
            args += ['--tracks', 'DNA']

        self.subprocess_check_call(args)

    def _add_json(self, json_data):
        cmd = [
            'perl', self._jbrowse_bin('add-json.pl'),
            json.dumps(json_data),
            os.path.join('data', 'trackList.json')
        ]
        self.subprocess_check_call(cmd)

    def _add_track_json(self, json_data):
        if len(json_data) == 0:
            return

        tmp = tempfile.NamedTemporaryFile(delete=False)
        json.dump(json_data, tmp)
        tmp.close()
        cmd = ['perl', self._jbrowse_bin('add-track-json.pl'), tmp.name,
               os.path.join('data', 'trackList.json')]
        self.subprocess_check_call(cmd)
        os.unlink(tmp.name)

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

        dest = os.path.join(self.outdir, 'data', 'raw', trackData['label'] + '.gff')

        self._sort_gff(gff3, dest)

        url = os.path.join('raw', trackData['label'] + '.gff.gz')
        trackData.update({
            "urlTemplate": url,
            "storeClass": "JBrowse/Store/SeqFeature/GFF3Tabix",
        })

        trackData['glyph'] = 'JBrowse/View/FeatureGlyph/Segments'

        trackData['trackType'] = 'BlastView/View/Track/CanvasFeatures'
        trackData['type'] = 'BlastView/View/Track/CanvasFeatures'

        self._add_track_json(trackData)

        os.unlink(gff3)

        if blastOpts.get('index', 'false') == 'true':
            self.tracksToIndex.append("%s" % trackData['label'])

    def add_bigwig(self, data, trackData, wiggleOpts, **kwargs):
        dest = os.path.join('data', 'raw', trackData['label'] + '.bw')
        self.symlink_or_copy(os.path.realpath(data), dest)

        url = os.path.join('raw', trackData['label'] + '.bw')
        trackData.update({
            "urlTemplate": url,
            "storeClass": "JBrowse/Store/SeqFeature/BigWig",
            "type": "JBrowse/View/Track/Wiggle/Density",
        })

        trackData['type'] = wiggleOpts['type']
        trackData['variance_band'] = True if wiggleOpts['variance_band'] == 'true' else False

        if 'min' in wiggleOpts and 'max' in wiggleOpts:
            trackData['min_score'] = wiggleOpts['min']
            trackData['max_score'] = wiggleOpts['max']
        else:
            trackData['autoscale'] = wiggleOpts.get('autoscale', 'local')

        trackData['scale'] = wiggleOpts['scale']

        self._add_track_json(trackData)

    def add_bigwig_multiple(self, data, trackData, wiggleOpts, **kwargs):

        urls = []
        for idx, bw in enumerate(data):
            dest = os.path.join('data', 'raw', trackData['label'] + '_' + str(idx) + '.bw')
            self.symlink_or_copy(bw[1], dest)

            urls.append({"url": os.path.join('raw', trackData['label'] + '_' + str(idx) + '.bw'), "name": str(idx + 1) + ' - ' + bw[0]})

        trackData.update({
            "urlTemplates": urls,
            "showTooltips": "true",
            "storeClass": "MultiBigWig/Store/SeqFeature/MultiBigWig",
            "type": "MultiBigWig/View/Track/MultiWiggle/MultiDensity",
        })
        if 'XYPlot' in wiggleOpts['type']:
            trackData['type'] = "MultiBigWig/View/Track/MultiWiggle/MultiXYPlot"

        trackData['variance_band'] = True if wiggleOpts['variance_band'] == 'true' else False

        if 'min' in wiggleOpts and 'max' in wiggleOpts:
            trackData['min_score'] = wiggleOpts['min']
            trackData['max_score'] = wiggleOpts['max']
        else:
            trackData['autoscale'] = wiggleOpts.get('autoscale', 'local')

        trackData['scale'] = wiggleOpts['scale']

        self._add_track_json(trackData)

    def add_maf(self, data, trackData, mafOpts, **kwargs):
        script = os.path.realpath(os.path.join(self.jbrowse, 'plugins', 'MAFViewer', 'bin', 'maf2bed.pl'))
        dest = os.path.join('data', 'raw', trackData['label'] + '.txt')

        tmp1 = tempfile.NamedTemporaryFile(delete=False)
        tmp1.close()

        # Process MAF to bed-like
        cmd = [script, data]
        self.subprocess_check_call(cmd, output=tmp1.path)

        # Sort / Index it
        self._sort_bed(tmp1.path, dest)
        # Cleanup
        try:
            os.remove(tmp1.path)
        except OSError:
            pass

        # Construct samples list
        # We could get this from galaxy metadata, not sure how easily.
        ps = subprocess.Popen(['grep', '^s [^ ]*', '-o', data], stdout=subprocess.PIPE)
        output = subprocess.check_output(('sort', '-u'), stdin=ps.stdout)
        ps.wait()
        samples = [x[2:] for x in output]

        trackData.update({
            "storeClass": "MAFViewer/Store/SeqFeature/MAFTabix",
            "type": "MAFViewer/View/Track/MAF",
            "urlTemplate": trackData['label'] + '.txt.gz',
            "samples": samples,
        })

        self._add_track_json(trackData)

    def add_bam(self, data, trackData, bamOpts, bam_index=None, **kwargs):
        dest = os.path.join('data', 'raw', trackData['label'] + '.bam')
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

        url = os.path.join('raw', trackData['label'] + '.bam')
        trackData.update({
            "urlTemplate": url,
            "type": "JBrowse/View/Track/Alignments2",
            "storeClass": "JBrowse/Store/SeqFeature/BAM",
            "chunkSizeLimit": bamOpts.get('chunkSizeLimit', '5000000')
        })

        # Apollo will only switch to the (prettier) 'bam-read' className if it's not set explicitly in the track config
        # So remove the default 'feature' value for these bam tracks
        if 'className' in trackData['style'] and trackData['style']['className'] == 'feature':
            del trackData['style']['className']

        self._add_track_json(trackData)

        if bamOpts.get('auto_snp', 'false') == 'true':
            trackData2 = copy.copy(trackData)
            trackData2.update({
                "type": "JBrowse/View/Track/SNPCoverage",
                "key": trackData['key'] + " - SNPs/Coverage",
                "label": trackData['label'] + "_autosnp",
                "chunkSizeLimit": bamOpts.get('chunkSizeLimit', '5000000')
            })
            self._add_track_json(trackData2)

    def add_vcf(self, data, trackData, vcfOpts={}, **kwargs):
        dest = os.path.join('data', 'raw', trackData['label'] + '.vcf')
        # ln?
        cmd = ['ln', '-s', data, dest]
        self.subprocess_check_call(cmd)
        cmd = ['bgzip', dest]
        self.subprocess_check_call(cmd)
        cmd = ['tabix', '-p', 'vcf', dest + '.gz']
        self.subprocess_check_call(cmd)

        url = os.path.join('raw', trackData['label'] + '.vcf.gz')
        trackData.update({
            "urlTemplate": url,
            "type": "JBrowse/View/Track/HTMLVariants",
            "storeClass": "JBrowse/Store/SeqFeature/VCFTabix",
        })
        self._add_track_json(trackData)

    def _sort_gff(self, data, dest):
        # Only index if not already done
        if not os.path.exists(dest):
            cmd = "gff3sort.pl --precise '%s' | grep -v \"^$\" > '%s'" % (data, dest)
            self.subprocess_popen(cmd)

            self.subprocess_check_call(['bgzip', '-f', dest])
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
        dest = os.path.join(self.outdir, 'data', 'raw', trackData['label'] + '.gff')

        self._sort_gff(data, dest)

        url = os.path.join('raw', trackData['label'] + '.gff.gz')
        trackData.update({
            "urlTemplate": url,
            "storeClass": "JBrowse/Store/SeqFeature/GFF3Tabix",
        })

        if 'match' in gffOpts:
            trackData['glyph'] = 'JBrowse/View/FeatureGlyph/Segments'

        trackType = 'JBrowse/View/Track/CanvasFeatures'
        if 'trackType' in gffOpts:
            trackType = gffOpts['trackType']
        trackData['type'] = trackType
        trackData['trackType'] = trackType  # Probably only used by old jbrowse versions

        if trackType in ['JBrowse/View/Track/CanvasFeatures', 'NeatCanvasFeatures/View/Track/NeatFeatures']:
            if 'transcriptType' in gffOpts and gffOpts['transcriptType']:
                trackData['transcriptType'] = gffOpts['transcriptType']
            if 'subParts' in gffOpts and gffOpts['subParts']:
                trackData['subParts'] = gffOpts['subParts']
            if 'impliedUTRs' in gffOpts and gffOpts['impliedUTRs']:
                trackData['impliedUTRs'] = gffOpts['impliedUTRs']
        elif trackType in ['JBrowse/View/Track/HTMLFeatures', 'NeatHTMLFeatures/View/Track/NeatFeatures']:
            if 'topLevelFeatures' in gffOpts and gffOpts['topLevelFeatures']:
                trackData['topLevelFeatures'] = gffOpts['topLevelFeatures']

        self._add_track_json(trackData)

        if gffOpts.get('index', 'false') == 'true':
            self.tracksToIndex.append("%s" % trackData['label'])

    def add_bed(self, data, format, trackData, gffOpts, **kwargs):
        dest = os.path.join(self.outdir, 'data', 'raw', trackData['label'] + '.bed')

        self._sort_bed(data, dest)

        url = os.path.join('raw', trackData['label'] + '.bed.gz')
        trackData.update({
            "urlTemplate": url,
            "storeClass": "JBrowse/Store/SeqFeature/BEDTabix",
        })

        if 'match' in gffOpts:
            trackData['glyph'] = 'JBrowse/View/FeatureGlyph/Segments'

        trackType = gffOpts.get('trackType', 'JBrowse/View/Track/CanvasFeatures')
        trackData['type'] = trackType

        if trackType in ['JBrowse/View/Track/CanvasFeatures', 'NeatCanvasFeatures/View/Track/NeatFeatures']:
            if 'transcriptType' in gffOpts and gffOpts['transcriptType']:
                trackData['transcriptType'] = gffOpts['transcriptType']
            if 'subParts' in gffOpts and gffOpts['subParts']:
                trackData['subParts'] = gffOpts['subParts']
            if 'impliedUTRs' in gffOpts and gffOpts['impliedUTRs']:
                trackData['impliedUTRs'] = gffOpts['impliedUTRs']
        elif trackType in ['JBrowse/View/Track/HTMLFeatures', 'NeatHTMLFeatures/View/Track/NeatFeatures']:
            if 'topLevelFeatures' in gffOpts and gffOpts['topLevelFeatures']:
                trackData['topLevelFeatures'] = gffOpts['topLevelFeatures']

        self._add_track_json(trackData)

        if gffOpts.get('index', 'false') == 'true':
            self.tracksToIndex.append("%s" % trackData['label'])

    def add_genbank(self, data, format, trackData, gffOpts, **kwargs):
        cmd = [
            'perl', self._jbrowse_bin('flatfile-to-json.pl'),
            '--genbank', data,
            '--trackLabel', trackData['label'],
            '--key', trackData['key']
        ]

        # className in --clientConfig is ignored, it needs to be set with --className
        if 'className' in trackData['style']:
            cmd += ['--className', trackData['style']['className']]

        config = copy.copy(trackData)
        clientConfig = trackData['style']
        del config['style']

        if 'match' in gffOpts:
            config['glyph'] = 'JBrowse/View/FeatureGlyph/Segments'
            if bool(gffOpts['match']):
                # Can be empty for CanvasFeatures = will take all by default
                cmd += ['--type', gffOpts['match']]

        cmd += ['--clientConfig', json.dumps(clientConfig)]

        trackType = 'JBrowse/View/Track/CanvasFeatures'
        if 'trackType' in gffOpts:
            trackType = gffOpts['trackType']

        if trackType == 'JBrowse/View/Track/CanvasFeatures':
            if 'transcriptType' in gffOpts and gffOpts['transcriptType']:
                config['transcriptType'] = gffOpts['transcriptType']
            if 'subParts' in gffOpts and gffOpts['subParts']:
                config['subParts'] = gffOpts['subParts']
            if 'impliedUTRs' in gffOpts and gffOpts['impliedUTRs']:
                config['impliedUTRs'] = gffOpts['impliedUTRs']
        elif trackType == 'JBrowse/View/Track/HTMLFeatures':
            if 'transcriptType' in gffOpts and gffOpts['transcriptType']:
                cmd += ['--type', gffOpts['transcriptType']]

        cmd += [
            '--trackType', gffOpts['trackType']
        ]

        cmd.extend(['--config', json.dumps(config)])

        self.subprocess_check_call(cmd)

        if gffOpts.get('index', 'false') == 'true':
            self.tracksToIndex.append("%s" % trackData['label'])

    def add_rest(self, url, trackData):
        data = {
            "label": trackData['label'],
            "key": trackData['key'],
            "category": trackData['category'],
            "type": "JBrowse/View/Track/HTMLFeatures",
            "storeClass": "JBrowse/Store/SeqFeature/REST",
            "baseUrl": url
        }
        self._add_track_json(data)

    def add_sparql(self, url, query, trackData):
        data = {
            "label": trackData['label'],
            "key": trackData['key'],
            "category": trackData['category'],
            "type": "JBrowse/View/Track/CanvasFeatures",
            "storeClass": "JBrowse/Store/SeqFeature/SPARQL",
            "urlTemplate": url,
            "queryTemplate": query
        }
        self._add_track_json(data)

    def traverse_to_option_parent(self, splitKey, outputTrackConfig):
        trackConfigSubDict = outputTrackConfig
        for part in splitKey[:-1]:
            if trackConfigSubDict.get(part) is None:
                trackConfigSubDict[part] = dict()
            trackConfigSubDict = trackConfigSubDict[part]
        assert isinstance(trackConfigSubDict, dict), 'Config element {} is not a dict'.format(trackConfigSubDict)
        return trackConfigSubDict

    def get_formatted_option(self, valType2ValDict, mapped_chars):
        assert isinstance(valType2ValDict, dict) and len(valType2ValDict.items()) == 1
        for valType, value in valType2ValDict.items():
            if valType == "text":
                for char, mapped_char in mapped_chars.items():
                    value = value.replace(mapped_char, char)
            elif valType == "integer":
                value = int(value)
            elif valType == "float":
                value = float(value)
            else:  # boolean
                value = {'true': True, 'false': False}[value]
            return value

    def set_custom_track_options(self, customTrackConfig, outputTrackConfig, mapped_chars):
        for optKey, optType2ValDict in customTrackConfig.items():
            splitKey = optKey.split('.')
            trackConfigOptionParent = self.traverse_to_option_parent(splitKey, outputTrackConfig)
            optVal = self.get_formatted_option(optType2ValDict, mapped_chars)
            trackConfigOptionParent[splitKey[-1]] = optVal

    def process_annotations(self, track):
        category = track['category'].replace('__pd__date__pd__', TODAY)
        outputTrackConfig = {
            'style': {
                'label': track['style'].get('label', 'description'),
                'className': track['style'].get('className', 'feature'),
                'description': track['style'].get('description', ''),
            },
            'overridePlugins': track['style'].get('overridePlugins', False) == 'True',
            'overrideDraggable': track['style'].get('overrideDraggable', False) == 'True',
            'maxHeight': track['style'].get('maxHeight', '600'),
            'category': category,
        }

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

            # Colour parsing is complex due to different track types having
            # different colour options.
            colourOptions = self.cs.parse_colours(track['conf']['options'], track['format'], gff3=dataset_path)
            # This used to be done with a dict.update() call, however that wiped out any previous style settings...
            for key in colourOptions:
                if key == 'style':
                    for subkey in colourOptions['style']:
                        outputTrackConfig['style'][subkey] = colourOptions['style'][subkey]
                else:
                    outputTrackConfig[key] = colourOptions[key]

            if 'menus' in track['conf']['options']:
                menus = self.cs.parse_menus(track['conf']['options'])
                outputTrackConfig.update(menus)

            customTrackConfig = track['conf']['options'].get('custom_config', {})
            if customTrackConfig:
                self.set_custom_track_options(customTrackConfig, outputTrackConfig, mapped_chars)

            # import pprint; pprint.pprint(track)
            # import sys; sys.exit()
            if dataset_ext in ('gff', 'gff3'):
                self.add_gff(dataset_path, dataset_ext, outputTrackConfig,
                             track['conf']['options']['gff'])
            elif dataset_ext in ('bed', ):
                self.add_bed(dataset_path, dataset_ext, outputTrackConfig,
                             track['conf']['options']['gff'])
            elif dataset_ext in ('genbank', ):
                self.add_genbank(dataset_path, dataset_ext, outputTrackConfig,
                                 track['conf']['options']['gff'])
            elif dataset_ext == 'bigwig':
                self.add_bigwig(dataset_path, outputTrackConfig,
                                track['conf']['options']['wiggle'])
            elif dataset_ext == 'bigwig_multiple':
                self.add_bigwig_multiple(dataset_path, outputTrackConfig,
                                         track['conf']['options']['wiggle'])
            elif dataset_ext == 'maf':
                self.add_maf(dataset_path, outputTrackConfig,
                             track['conf']['options']['maf'])
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
            elif dataset_ext == 'rest':
                self.add_rest(track['conf']['options']['rest']['url'], outputTrackConfig)
            elif dataset_ext == 'sparql':
                sparql_query = track['conf']['options']['sparql']['query']
                for key, value in mapped_chars.items():
                    sparql_query = sparql_query.replace(value, key)
                self.add_sparql(track['conf']['options']['sparql']['url'], sparql_query, outputTrackConfig)
            else:
                log.warn('Do not know how to handle %s', dataset_ext)

            # Return non-human label for use in other fields
            yield outputTrackConfig['label']

    def add_final_data(self, data):
        viz_data = {}
        if len(data['visibility']['default_on']) > 0:
            viz_data['defaultTracks'] = ','.join(data['visibility']['default_on'])

        if len(data['visibility']['always']) > 0:
            viz_data['alwaysOnTracks'] = ','.join(data['visibility']['always'])

        if len(data['visibility']['force']) > 0:
            viz_data['forceTracks'] = ','.join(data['visibility']['force'])

        generalData = {}
        if data['general']['aboutDescription'] is not None:
            generalData['aboutThisBrowser'] = {'description': data['general']['aboutDescription'].strip()}

        generalData['view'] = {
            'trackPadding': data['general']['trackPadding']
        }
        generalData['shareLink'] = (data['general']['shareLink'] == 'true')
        generalData['show_tracklist'] = (data['general']['show_tracklist'] == 'true')
        generalData['show_nav'] = (data['general']['show_nav'] == 'true')
        generalData['show_overview'] = (data['general']['show_overview'] == 'true')
        generalData['show_menu'] = (data['general']['show_menu'] == 'true')
        generalData['hideGenomeOptions'] = (data['general']['hideGenomeOptions'] == 'true')
        generalData['plugins'] = data['plugins']

        viz_data.update(generalData)
        self._add_json(viz_data)

        if 'GCContent' in data['plugins_python']:
            self._add_track_json({
                "storeClass": "JBrowse/Store/Sequence/IndexedFasta",
                "type": "GCContent/View/Track/GCContentXY",
                "label": "GC Content",
                "key": "GCContentXY",
                "urlTemplate": "seq/genome.fasta",
                "bicolor_pivot": 0.5,
                "category": "GC Content",
                "metadata": {
                    "tool_tool": '<a target="_blank" href="https://github.com/elsiklab/gccontent/commit/030180e75a19fad79478d43a67c566ec6">elsiklab/gccontent</a>',
                    "tool_tool_version": "5c8b0582ecebf9edf684c76af8075fb3d30ec3fa",
                    "dataset_edam_format": "",
                    "dataset_size": "",
                    "history_display_name": "",
                    "history_user_email": "",
                    "metadata_dbkey": "",
                }
                # TODO: Expose params for everyone.
            })
            self._add_track_json({
                "storeClass": "JBrowse/Store/Sequence/IndexedFasta",
                "type": "GCContent/View/Track/GCContentXY",
                "label": "GC skew",
                "key": "GCSkew",
                "urlTemplate": "seq/genome.fasta",
                "gcMode": "skew",
                "min_score": -1,
                "bicolor_pivot": 0,
                "category": "GC Content",
                "metadata": {
                    "tool_tool": '<a target="_blank" href="https://github.com/elsiklab/gccontent/commit/030180e75a19fad79478d43a67c566ec6">elsiklab/gccontent</a>',
                    "tool_tool_version": "5c8b0582ecebf9edf684c76af8075fb3d30ec3fa",
                    "dataset_edam_format": "",
                    "dataset_size": "",
                    "history_display_name": "",
                    "history_user_email": "",
                    "metadata_dbkey": "",
                }
                # TODO: Expose params for everyone.
            })

        if 'ComboTrackSelector' in data['plugins_python']:
            with open(os.path.join(self.outdir, 'data', 'trackList.json'), 'r') as handle:
                trackListJson = json.load(handle)
                trackListJson.update({
                    "trackSelector": {
                        "renameFacets": {
                            "tool_tool": "Tool ID",
                            "tool_tool_id": "Tool ID",
                            "tool_tool_version": "Tool Version",
                            "dataset_edam_format": "EDAM",
                            "dataset_size": "Size",
                            "history_display_name": "History Name",
                            "history_user_email": "Owner",
                            "metadata_dbkey": "Dbkey",
                        },
                        "displayColumns": [
                            "key",
                            "tool_tool",
                            "tool_tool_version",
                            "dataset_edam_format",
                            "dataset_size",
                            "history_display_name",
                            "history_user_email",
                            "metadata_dbkey",
                        ],
                        "type": "Faceted",
                        "title": ["Galaxy Metadata"],
                        "icon": "https://galaxyproject.org/images/logos/galaxy-icon-square.png",
                        "escapeHTMLInData": False
                    },
                    "trackMetadata": {
                        "indexFacets": [
                            "category",
                            "key",
                            "tool_tool_id",
                            "tool_tool_version",
                            "dataset_edam_format",
                            "history_user_email",
                            "history_display_name"
                        ]
                    }
                })
                with open(os.path.join(self.outdir, 'data', 'trackList2.json'), 'w') as handle:
                    json.dump(trackListJson, handle)

    def clone_jbrowse(self, jbrowse_dir, destination, minimal=False):
        """Clone a JBrowse directory into a destination directory.
        """
        if minimal:
            # Should be the absolute minimum required for JBrowse to function.
            interesting = [
                'dist', 'img', 'index.html', 'jbrowse.conf', 'jbrowse_conf.json', 'webpack.config.js'
            ]
            for i in interesting:
                cmd = ['cp', '-r', os.path.join(jbrowse_dir, i), destination]
                self.subprocess_check_call(cmd)
        else:
            # JBrowse seems to have included some bad symlinks, cp ignores bad symlinks
            # unlike copytree
            cmd = ['cp', '-r', os.path.join(jbrowse_dir, '.'), destination]
            self.subprocess_check_call(cmd)

        cmd = ['mkdir', '-p', os.path.join(destination, 'data', 'raw')]
        self.subprocess_check_call(cmd)

        # http://unix.stackexchange.com/a/38691/22785
        # JBrowse releases come with some broken symlinks
        cmd = ['find', destination, '-type', 'l', '-xtype', 'l']
        symlinks = self.subprocess_check_output(cmd)

        for i in symlinks:
            try:
                os.unlink(i)
            except OSError:
                pass


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
        'plugins': [],
        'plugins_python': [],
    }

    plugins = root.find('plugins').attrib
    if plugins['GCContent'] == 'True':
        extra_data['plugins_python'].append('GCContent')
        extra_data['plugins'].append({
            'location': 'https://cdn.jsdelivr.net/gh/elsiklab/gccontent@5c8b0582ecebf9edf684c76af8075fb3d30ec3fa/',
            'name': 'GCContent'
        })

    # Not needed in 1.16.1: it's built in the conda package now, and this plugin doesn't need to be enabled anywhere
    # if plugins['Bookmarks'] == 'True':
    #    extra_data['plugins'].append({
    #        'location': 'https://cdn.jsdelivr.net/gh/TAMU-CPT/bookmarks-jbrowse@5242694120274c86e1ccd5cb0e5e943e78f82393/',
    #        'name': 'Bookmarks'
    #    })

    # Not needed in 1.16.1: it's built in the conda package now, and this plugin doesn't need to be enabled anywhere
    if plugins['ComboTrackSelector'] == 'True':
        extra_data['plugins_python'].append('ComboTrackSelector')
        # Not needed in 1.16.1: it's built in the conda package now, and this plugin doesn't need to be enabled anywhere
        # extra_data['plugins'].append({
        #    'location': 'https://cdn.jsdelivr.net/gh/Arabidopsis-Information-Portal/ComboTrackSelector@52403928d5ccbe2e3a86b0fa5eb8e61c0f2e2f57/',
        #    'icon': 'https://galaxyproject.org/images/logos/galaxy-icon-square.png',
        #    'name': 'ComboTrackSelector'
        # })

    if plugins['theme'] == 'Minimalist':
        extra_data['plugins'].append({
            'location': 'https://cdn.jsdelivr.net/gh/erasche/jbrowse-minimalist-theme@d698718442da306cf87f033c72ddb745f3077775/',
            'name': 'MinimalistTheme'
        })
    elif plugins['theme'] == 'Dark':
        extra_data['plugins'].append({
            'location': 'https://cdn.jsdelivr.net/gh/erasche/jbrowse-dark-theme@689eceb7e33bbc1b9b15518d45a5a79b2e5d0a26/',
            'name': 'DarkTheme'
        })

    if plugins['BlastView'] == 'True':
        extra_data['plugins_python'].append('BlastView')
        extra_data['plugins'].append({
            'location': 'https://cdn.jsdelivr.net/gh/TAMU-CPT/blastview@97572a21b7f011c2b4d9a0b5af40e292d694cbef/',
            'name': 'BlastView'
        })

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
        keys = jc.process_annotations(track_conf)

        for key in keys:
            extra_data['visibility'][track.attrib.get('visibility', 'default_off')].append(key)

    jc.add_final_data(extra_data)
    jc.generate_names()
