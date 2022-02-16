#!/usr/bin/env python
import argparse
import binascii
import datetime
import hashlib
import logging
import os
import shutil
import struct
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from collections import defaultdict

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

    def __init__(self, jbrowse, outdir, genomes):
        self.cs = ColorScaling()
        self.jbrowse = jbrowse
        self.outdir = outdir
        self.genome_paths = genomes
        self.tracksToIndex = []

        self.clone_jbrowse(self.jbrowse, self.outdir)

        self.process_genomes()

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

    def symlink_or_copy(self, src, dest):
        if 'GALAXY_JBROWSE_SYMLINKS' in os.environ and bool(os.environ['GALAXY_JBROWSE_SYMLINKS']):
            cmd = ['ln', '-s', src, dest]
        else:
            cmd = ['cp', src, dest]

        return self.subprocess_check_call(cmd)

    def symlink_or_copy_load_action(self):
        if 'GALAXY_JBROWSE_SYMLINKS' in os.environ and bool(os.environ['GALAXY_JBROWSE_SYMLINKS']):
            return 'symlink'
        else:
            return 'copy'

    def process_genomes(self):
        for genome_node in self.genome_paths:
            # We only expect one input genome per run. This for loop is just
            # easier to write than the alternative / catches any possible
            # issues.
            self.add_assembly(genome_node['path'], genome_node['label'])

    def add_assembly(self, path, label):
        copied_genome = os.path.join(self.outdir, 'data', 'genome.fasta')
        shutil.copy(path, copied_genome)

        # Compress with bgzip
        cmd = ['bgzip', copied_genome]
        self.subprocess_check_call(cmd)

        # FAI Index
        cmd = ['samtools', 'faidx', copied_genome + '.gz']
        self.subprocess_check_call(cmd)

        self.subprocess_check_call([
            'jbrowse', 'add-assembly',
            '--load', 'inPlace',
            '--name', label,
            '--type', 'bgzipFasta',
            '--target', os.path.join(self.outdir, 'data'),
            '--skipCheck',
            os.path.join('data', 'genome.fasta.gz')])

    def text_index(self):
        # Index tracks
        args = [
            'jbrowse', 'text-index',
            '--target', os.path.join(self.outdir, 'data')
        ]

        tracks = ','.join(self.tracksToIndex)
        if tracks:
            args += ['--tracks', tracks]

            self.subprocess_check_call(args)

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

        rel_dest = os.path.join('data', trackData['label'] + '.gff')
        dest = os.path.join(self.outdir, rel_dest)

        self._sort_gff(gff3, dest)
        os.unlink(gff3)

        self._add_track(trackData['label'], trackData['key'], trackData['category'], rel_dest + '.gz')

    def add_bigwig(self, data, trackData, wiggleOpts, **kwargs):

        rel_dest = os.path.join('data', trackData['label'] + '.bw')
        dest = os.path.join(self.outdir, rel_dest)
        self.symlink_or_copy(os.path.realpath(data), dest)

        self._add_track(trackData['label'], trackData['key'], trackData['category'], rel_dest)

    # Anything ending in "am" (Bam or Cram)
    def add_xam(self, data, trackData, xamOpts, index=None, ext="bam", **kwargs):

        index_ext = "bai"
        if ext == "cram":
            index_ext = "crai"

        rel_dest = os.path.join('data', trackData['label'] + '.%s' % ext)
        dest = os.path.join(self.outdir, rel_dest)

        self.symlink_or_copy(os.path.realpath(data), dest)

        if index is not None and os.path.exists(os.path.realpath(index)):
            # xai most probably made by galaxy and stored in galaxy dirs, need to copy it to dest
            self.subprocess_check_call(['cp', os.path.realpath(index), dest + '.%s' % index_ext])
        else:
            # Can happen in exotic condition
            # e.g. if bam imported as symlink with datatype=unsorted.bam, then datatype changed to bam
            #      => no index generated by galaxy, but there might be one next to the symlink target
            #      this trick allows to skip the bam sorting made by galaxy if already done outside
            if os.path.exists(os.path.realpath(data) + '.%s' % index_ext):
                self.symlink_or_copy(os.path.realpath(data) + '.%s' % index_ext, dest + '.%s' % index_ext)
            else:
                log.warn('Could not find a bam index (.%s file) for %s', (index_ext, data))

        self._add_track(trackData['label'], trackData['key'], trackData['category'], rel_dest)

    def add_vcf(self, data, trackData, vcfOpts={}, zipped=False, **kwargs):

        if zipped:
            rel_dest = os.path.join('data', trackData['label'] + '.vcf.gz')
            dest = os.path.join(self.outdir, rel_dest)
            shutil.copy(os.path.realpath(data), dest)
        else:
            rel_dest = os.path.join('data', trackData['label'] + '.vcf')
            dest = os.path.join(self.outdir, rel_dest)
            shutil.copy(os.path.realpath(data), dest)

            cmd = ['bgzip', dest]
            self.subprocess_check_call(cmd)
            cmd = ['tabix', dest + '.gz']
            self.subprocess_check_call(cmd)

            rel_dest = os.path.join('data', trackData['label'] + '.vcf.gz')

        self._add_track(trackData['label'], trackData['key'], trackData['category'], rel_dest)

    def add_gff(self, data, format, trackData, gffOpts, **kwargs):
        rel_dest = os.path.join('data', trackData['label'] + '.gff')
        dest = os.path.join(self.outdir, rel_dest)

        self._sort_gff(data, dest)

        self._add_track(trackData['label'], trackData['key'], trackData['category'], rel_dest + '.gz')

    def add_bed(self, data, format, trackData, gffOpts, **kwargs):
        rel_dest = os.path.join('data', trackData['label'] + '.bed')
        dest = os.path.join(self.outdir, rel_dest)

        self._sort_bed(data, dest)

        self._add_track(trackData['label'], trackData['key'], trackData['category'], rel_dest + '.gz')

    def add_paf(self, data, trackData, pafOpts, **kwargs):
        rel_dest = os.path.join('data', trackData['label'] + '.paf')
        dest = os.path.join(self.outdir, rel_dest)

        self.symlink_or_copy(os.path.realpath(data), dest)

        self.add_assembly(pafOpts['genome'], pafOpts['genome_label'])

        self._add_track(trackData['label'], trackData['key'], trackData['category'], rel_dest)

    def add_sparql(self, url, query, query_refnames, trackData):

        json_track_data = {
            "type": "FeatureTrack",
            "trackId": id,
            "name": trackData['label'],
            "adapter": {
                "type": "SPARQLAdapter",
                "endpoint": {
                    "uri": url,
                    "locationType": "UriLocation"
                },
                "queryTemplate": query
            },
            "category": [
                trackData['category']
            ]
        }

        if query_refnames:
            json_track_data['adapter']['refNamesQueryTemplate']: query_refnames

        self.subprocess_check_call([
            'jbrowse', 'add-track-json',
            '--target', os.path.join(self.outdir, 'data'),
            json_track_data])

        # Doesn't work as of 1.6.4, might work in the future
        # self.subprocess_check_call([
        #     'jbrowse', 'add-track',
        #     '--trackType', 'sparql',
        #     '--name', trackData['label'],
        #     '--category', trackData['category'],
        #     '--target', os.path.join(self.outdir, 'data'),
        #     '--trackId', id,
        #     '--config', '{"queryTemplate": "%s"}' % query,
        #     url])

    def _add_track(self, id, label, category, path):
        self.subprocess_check_call([
            'jbrowse', 'add-track',
            '--load', 'inPlace',
            '--name', label,
            '--category', category,
            '--target', os.path.join(self.outdir, 'data'),
            '--trackId', id,
            path])

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
                'description': track['style'].get('description', ''),
            },
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
            elif dataset_ext == 'bed':
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

                self.add_xam(dataset_path, outputTrackConfig,
                             track['conf']['options']['pileup'],
                             index=real_indexes[i], ext="bam")
            elif dataset_ext == 'cram':
                real_indexes = track['conf']['options']['cram']['cram_indices']['cram_index']
                if not isinstance(real_indexes, list):
                    # <bam_indices>
                    #  <bam_index>/path/to/a.bam.bai</bam_index>
                    # </bam_indices>
                    #
                    # The above will result in the 'bam_index' key containing a
                    # string. If there are two or more indices, the container
                    # becomes a list. Fun!
                    real_indexes = [real_indexes]

                self.add_xam(dataset_path, outputTrackConfig,
                             track['conf']['options']['cram'],
                             index=real_indexes[i], ext="cram")
            elif dataset_ext == 'blastxml':
                self.add_blastxml(dataset_path, outputTrackConfig, track['conf']['options']['blast'])
            elif dataset_ext == 'vcf':
                self.add_vcf(dataset_path, outputTrackConfig)
            elif dataset_ext == 'vcf_bgzip':
                self.add_vcf(dataset_path, outputTrackConfig, zipped=True)
            elif dataset_ext == 'rest':
                self.add_rest(track['conf']['options']['rest']['url'], outputTrackConfig)
            elif dataset_ext == 'paf':
                self.add_paf(dataset_path, outputTrackConfig,
                             track['conf']['options']['paf'])
            elif dataset_ext == 'sparql':
                sparql_query = track['conf']['options']['sparql']['query']
                for key, value in mapped_chars.items():
                    sparql_query = sparql_query.replace(value, key)
                sparql_query_refnames = track['conf']['options']['sparql']['query_refnames']
                for key, value in mapped_chars.items():
                    sparql_query_refnames = sparql_query_refnames.replace(value, key)
                self.add_sparql(track['conf']['options']['sparql']['url'], sparql_query, sparql_query_refnames, outputTrackConfig)
            else:
                log.warn('Do not know how to handle %s', dataset_ext)

            # Return non-human label for use in other fields
            yield outputTrackConfig['label']

    def clone_jbrowse(self, jbrowse_dir, destination):
        """Clone a JBrowse directory into a destination directory.
        """

        copytree(jbrowse_dir, destination)

        try:
            shutil.rmtree(os.path.join(destination, 'test_data'))
        except OSError as e:
            print("Error: %s - %s." % (e.filename, e.strerror))

        os.makedirs(os.path.join(destination, 'data'))
        print("makedir %s" % (os.path.join(destination, 'data')))

        os.symlink('./data/config.json', os.path.join(destination, 'config.json'))


def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument('xml', type=argparse.FileType('r'), help='Track Configuration')

    parser.add_argument('--jbrowse', help='Folder containing a jbrowse release')
    parser.add_argument('--outdir', help='Output directory', default='out')
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
                'meta': metadata_from_node(x.find('metadata')),
                'label': x.attrib['label']
            }
            for x in root.findall('metadata/genomes/genome')
        ]
    )

    # TODO is this still needed?
    for track in root.findall('tracks/track'):
        track_conf = {}
        track_conf['trackfiles'] = []

        trackfiles = track.findall('files/trackFile')
        if trackfiles:
            for x in track.findall('files/trackFile'):
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

    jc.text_index()
