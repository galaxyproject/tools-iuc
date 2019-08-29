#!/usr/bin/env python
import argparse
import copy
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


class ColorScaling(object):

    COLOR_FUNCTION_TEMPLATE = """
    function(feature, variableName, glyphObject, track) {{
        var score = {score};
        {opacity}
        return 'rgba({red}, {green}, {blue}, ' + opacity + ')';
    }}
    """

    COLOR_FUNCTION_TEMPLATE_QUAL = """
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
            var opacity = (score - ({min})) / (({max}) - ({min}));
            opacity = Math.log10(opacity) + Math.log10({max});
        """,
        'blast': """
            var opacity = 0;
            if(score == 0.0) {
                opacity = 1;
            } else{
                opacity = (20 - Math.log10(score)) / 180;
            }
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
        return struct.unpack('BBB', hexstr.decode('hex'))

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
        trackConfig = {'menuTemplate': [{}, {}, {}]}

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
    d = {t.tag: {} if t.attrib else None}
    children = list(t)
    if children:
        dd = defaultdict(list)
        for dc in map(etree_to_dict, children):
            for k, v in dc.iteritems():
                dd[k].append(v)
        d = {t.tag: {k: v[0] if len(v) == 1 else v for k, v in dd.iteritems()}}
    if t.attrib:
        d[t.tag].update(('@' + k, v) for k, v in t.attrib.iteritems())
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


class JbrowseConnector(object):

    def __init__(self, jbrowse, outdir, genomes, standalone=False, gencode=1):
        self.TN_TABLE = {
            'gff3': '--gff',
            'gff': '--gff',
            'bed': '--bed',
            'genbank': '--gbk',
        }

        self.cs = ColorScaling()
        self.jbrowse = jbrowse
        self.outdir = outdir
        self.genome_paths = genomes
        self.standalone = standalone
        self.gencode = gencode
        self.tracksToIndex = []

        if standalone:
            self.clone_jbrowse(self.jbrowse, self.outdir)
        else:
            try:
                os.makedirs(self.outdir)
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

    def subprocess_check_call(self, command):
        log.debug('cd %s && %s', self.outdir, ' '.join(command))
        subprocess.check_call(command, cwd=self.outdir)

    def _jbrowse_bin(self, command):
        return os.path.realpath(os.path.join(self.jbrowse, 'bin', command))

    def process_genomes(self):
        for genome_path in self.genome_paths:
            self.subprocess_check_call([
                'perl', self._jbrowse_bin('prepare-refseqs.pl'),
                '--fasta', genome_path])

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
        if len(json_data.keys()) == 0:
            return

        tmp = tempfile.NamedTemporaryFile(delete=False)
        tmp.write(json.dumps(json_data))
        tmp.close()
        cmd = ['perl', self._jbrowse_bin('add-track-json.pl'), tmp.name,
               os.path.join('data', 'trackList.json')]
        self.subprocess_check_call(cmd)
        os.unlink(tmp.name)

    def _blastxml_to_gff3(self, xml, min_gap=10):
        gff3_unrebased = tempfile.NamedTemporaryFile(delete=False)
        cmd = ['python', os.path.join(INSTALLED_TO, 'blastxml_to_gapped_gff3.py'),
               '--trim', '--trim_end', '--min_gap', str(min_gap), xml]
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

        config = {
            'glyph': 'JBrowse/View/FeatureGlyph/Segments',
            "category": trackData['category'],
        }

        clientConfig = trackData['style']

        cmd = ['perl', self._jbrowse_bin('flatfile-to-json.pl'),
               '--gff', gff3,
               '--trackLabel', trackData['label'],
               '--key', trackData['key'],
               '--clientConfig', json.dumps(clientConfig),
               '--config', json.dumps(config),
               '--trackType', 'JBrowse/View/Track/CanvasFeatures'
               ]

        self.subprocess_check_call(cmd)
        os.unlink(gff3)

        if blastOpts.get('index', 'false') == 'true':
            self.tracksToIndex.append("%s" % trackData['label'])

    def add_bigwig(self, data, trackData, wiggleOpts, **kwargs):
        dest = os.path.join('data', 'raw', trackData['label'] + '.bw')
        cmd = ['ln', '-s', data, dest]
        self.subprocess_check_call(cmd)

        trackData.update({
            "urlTemplate": os.path.join('..', dest),
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

        self._add_track_json(trackData)

    def add_bam(self, data, trackData, bamOpts, bam_index=None, **kwargs):
        dest = os.path.join('data', 'raw', trackData['label'] + '.bam')
        cmd = ['ln', '-s', os.path.realpath(data), dest]
        self.subprocess_check_call(cmd)

        cmd = ['ln', '-s', os.path.realpath(bam_index), dest + '.bai']
        self.subprocess_check_call(cmd)

        trackData.update({
            "urlTemplate": os.path.join('..', dest),
            "type": "JBrowse/View/Track/Alignments2",
            "storeClass": "JBrowse/Store/SeqFeature/BAM",
        })

        self._add_track_json(trackData)

        if bamOpts.get('auto_snp', 'false') == 'true':
            trackData2 = copy.copy(trackData)
            trackData2.update({
                "type": "JBrowse/View/Track/SNPCoverage",
                "key": trackData['key'] + " - SNPs/Coverage",
                "label": trackData['label'] + "_autosnp",
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

        trackData.update({
            "urlTemplate": os.path.join('..', dest + '.gz'),
            "type": "JBrowse/View/Track/HTMLVariants",
            "storeClass": "JBrowse/Store/SeqFeature/VCFTabix",
        })
        self._add_track_json(trackData)

    def add_features(self, data, data_format, trackData, gffOpts, **kwargs):
        cmd = [
            'perl', self._jbrowse_bin('flatfile-to-json.pl'),
            self.TN_TABLE[data_format],
            data,
            '--trackLabel', trackData['label'],
            # '--trackType', 'JBrowse/View/Track/CanvasFeatures',
            '--key', trackData['key']
        ]

        config = copy.copy(trackData)
        clientConfig = trackData['style']
        del config['style']

        if 'match' in gffOpts:
            config['glyph'] = 'JBrowse/View/FeatureGlyph/Segments'
            cmd += ['--type', gffOpts['match']]

        cmd += ['--clientConfig', json.dumps(clientConfig),
                ]

        if 'trackType' in gffOpts:
            cmd += [
                '--trackType', gffOpts['trackType']
            ]
        else:
            cmd += [
                '--trackType', 'JBrowse/View/Track/CanvasFeatures'
            ]

        cmd.extend(['--config', json.dumps(config)])

        self.subprocess_check_call(cmd)

        if gffOpts.get('index', 'false') == 'true':
            self.tracksToIndex.append("%s" % trackData['label'])

    def process_annotations(self, track):
        outputTrackConfig = {
            'style': {
                'label': track['style'].get('label', 'description'),
                'className': track['style'].get('className', 'feature'),
                'description': track['style'].get('description', ''),
            },
            'category': track['category'],
        }

        for i, (dataset_path, dataset_ext, track_human_label) in enumerate(track['trackfiles']):
            log.info('Processing %s / %s', track['category'], track_human_label)
            outputTrackConfig['key'] = track_human_label
            hashData = [dataset_path, track_human_label, track['category']]
            outputTrackConfig['label'] = hashlib.md5('|'.join(hashData)).hexdigest() + '_%s' % i

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

            menus = self.cs.parse_menus(track['conf']['options'])
            outputTrackConfig.update(menus)

            # import pprint; pprint.pprint(track)
            # import sys; sys.exit()
            if dataset_ext in ('gff', 'gff3', 'bed', 'genbank'):
                self.add_features(dataset_path, dataset_ext, outputTrackConfig,
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

        viz_data.update(generalData)
        self._add_json(viz_data)

    def clone_jbrowse(self, jbrowse_dir, destination):
        """Clone a JBrowse directory into a destination directory.
        """
        # JBrowse seems to have included some bad symlinks, cp ignores bad symlinks
        # unlike copytree
        cmd = ['cp', '-r', os.path.join(jbrowse_dir, '.'), destination]
        log.debug(' '.join(cmd))
        subprocess.check_call(cmd)
        cmd = ['mkdir', '-p', os.path.join(destination, 'data', 'raw')]
        log.debug(' '.join(cmd))
        subprocess.check_call(cmd)

        # http://unix.stackexchange.com/a/38691/22785
        # JBrowse releases come with some broken symlinks
        cmd = ['find', destination, '-type', 'l', '-xtype', 'l', '-exec', 'rm', '{}', '+']
        log.debug(' '.join(cmd))
        subprocess.check_call(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument("xml", type=argparse.FileType("r"), help="Track Configuration")

    parser.add_argument('--jbrowse', help='Folder containing a jbrowse release')
    parser.add_argument('--outdir', help='Output directory', default='out')
    parser.add_argument('--standalone', help='Standalone mode includes a copy of JBrowse', action='store_true')
    args = parser.parse_args()

    tree = ET.parse(args.xml.name)
    root = tree.getroot()

    jc = JbrowseConnector(
        jbrowse=args.jbrowse,
        outdir=args.outdir,
        genomes=[os.path.realpath(x.attrib['path']) for x in root.findall('metadata/genomes/genome')],
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
        }
    }
    for track in root.findall('tracks/track'):
        track_conf = {}
        track_conf['trackfiles'] = [
            (os.path.realpath(x.attrib['path']), x.attrib['ext'], x.attrib['label'])
            for x in track.findall('files/trackFile')
        ]

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
