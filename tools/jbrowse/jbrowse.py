#!/usr/bin/env python
from string import Template
import os
import argparse
import subprocess
import hashlib
import tempfile
import json
import yaml
import logging
logging.basicConfig()
log = logging.getLogger(__name__)

COLOR_FUNCTION_TEMPLATE = Template("""
function(feature, variableName, glyphObject, track) {
    var score = ${score};
    ${opacity}
    return 'rgba(${red}, ${green}, ${blue}, ' + opacity + ')';
}
""")

BLAST_OPACITY_MATH = """
var opacity = 0;
if(score == 0.0) {
    opacity = 1;
} else{
    opacity = (20 - Math.log10(score)) / 180;
}
"""

BREWER_COLOUR_IDX = 0
BREWER_COLOUR_SCHEMES = [
    (228, 26, 28),
    (55, 126, 184),
    (77, 175, 74),
    (152, 78, 163),
    (255, 127, 0),
]


# score comes from feature._parent.get('score') or feature.get('score')
# Opacity math

TN_TABLE = {
    'gff3': '--gff',
    'gff': '--gff',
    'bed': '--bed',
    # 'genbank': '--gbk',
}

INSTALLED_TO = os.path.dirname(os.path.realpath(__file__))


class JbrowseConnector(object):

    def __init__(self, jbrowse, jbrowse_dir, outdir, genome):
        self.jbrowse = jbrowse
        self.jbrowse_dir = jbrowse_dir
        self.outdir = outdir
        self.genome_path = genome
        self.brewer_colour_idx = 0

        self.clone_jbrowse(self.jbrowse, self.outdir)
        self.process_genome()

    def subprocess_check_call(self, command):
        log.debug('cd %s && %s', self.jbrowse_dir, ' '.join(command))
        subprocess.check_call(command, cwd=self.jbrowse_dir)

    def _get_colours(self):
        r, g, b = BREWER_COLOUR_SCHEMES[self.brewer_colour_idx]
        self.brewer_colour_idx += 1
        return r, g, b

    def process_genome(self):
        self.subprocess_check_call(['perl', 'bin/prepare-refseqs.pl',
                                    '--fasta', self.genome_path])

    def _add_json(self, json_data):
        if len(json_data.keys()) == 0:
            return

        tmp = tempfile.NamedTemporaryFile(delete=False)
        tmp.write(json.dumps(json_data))
        tmp.close()
        cmd = ['perl', 'bin/add-track-json.pl', tmp.name,
               os.path.join('data', 'trackList.json')]
        self.subprocess_check_call(cmd)
        os.unlink(tmp.name)

    def add_blastxml(self, data, key, **kwargs):
        gff3_unrebased = tempfile.NamedTemporaryFile(delete=False)
        cmd = ['python', os.path.join(INSTALLED_TO, 'blastxml_to_gapped_gff3.py'),
               '--trim_end', '--min_gap', str(kwargs['min_gap']), data]
        subprocess.check_call(cmd, cwd=self.jbrowse_dir, stdout=gff3_unrebased)
        gff3_unrebased.close()

        gff3_rebased = tempfile.NamedTemporaryFile(delete=False)
        cmd = ['python', os.path.join(INSTALLED_TO, 'gff3_rebase.py')]
        if kwargs['protein']:
            cmd.append('--protein2dna')
        cmd.extend([kwargs['parent'], gff3_unrebased.name])
        subprocess.check_call(cmd, cwd=self.jbrowse_dir, stdout=gff3_rebased)
        gff3_rebased.close()

        label = hashlib.md5(data).hexdigest()

        red, green, blue = self._get_colours()
        color_function = COLOR_FUNCTION_TEMPLATE.substitute({
            'score': "feature._parent.get('score')",
            'opacity': BLAST_OPACITY_MATH,
            'red': red,
            'green': green,
            'blue': blue,
        })

        clientConfig = {
            'label': 'description',
            'color': color_function.replace('\n', ''),
            'description': 'Hit_titles',
        }
        config = {'glyph': 'JBrowse/View/FeatureGlyph/Segments'}
        if 'category' in kwargs:
            config['category'] = kwargs['category']

        cmd = ['perl', 'bin/flatfile-to-json.pl',
               '--gff', gff3_rebased.name,
               '--trackLabel', label,
               '--key', key,
               '--clientConfig', json.dumps(clientConfig),
               '--config', json.dumps(config),
               '--trackType', 'JBrowse/View/Track/CanvasFeatures'
               ]

        self.subprocess_check_call(cmd)
        os.unlink(gff3_rebased.name)
        os.unlink(gff3_unrebased.name)

    def _min_max_gff(self, gff_file):
        min_val = None
        max_val = None
        with open(gff_file, 'r') as handle:
            for line in handle:
                try:
                    value = float(line.split('\t')[5])
                    if min_val is None:
                        min_val = value
                    if max_val is None:
                        max_val = value

                    if value < min_val:
                        min_val = value

                    if value > max_val:
                        max_val = value
                except Exception:
                    pass
        return min_val, max_val

    def add_bigwig(self, data, key, **kwargs):
        label = hashlib.md5(data).hexdigest()
        dest = os.path.join('data', 'raw', os.path.basename(data))
        # ln?
        cmd = ['ln', '-s', data, dest]
        self.subprocess_check_call(cmd)

        track_data = {
            "label": label,
            "urlTemplate": os.path.join('..', dest),
            "bicolor_pivot": "zero",
            "storeClass": "JBrowse/Store/SeqFeature/BigWig",
            "type": "JBrowse/View/Track/Wiggle/Density",
            "key": key,
        }
        track_data.update(kwargs)

        if 'min' not in track_data and 'max' not in track_data \
                and 'autoscale' not in track_data:
            track_data['autoscale'] = 'local'

        self._add_json(track_data)

    def add_bam(self, data, key, **kwargs):
        label = hashlib.md5(data).hexdigest()
        dest = os.path.join('data', 'raw', os.path.basename(data))
        # ln?
        cmd = ['ln', '-s', data, dest]
        self.subprocess_check_call(cmd)

        bai_source = kwargs['bam_index']
        cmd = ['ln', '-s', bai_source, dest + '.bai']
        self.subprocess_check_call(cmd)

        track_data = {
            "urlTemplate": os.path.join('..', dest),
            "key": key,
            "label": label,
            "type": "JBrowse/View/Track/Alignments2",
            "storeClass": "JBrowse/Store/SeqFeature/BAM",
        }
        if 'category' in kwargs:
            track_data['category'] = kwargs['category']
        self._add_json(track_data)

        if kwargs.get('auto_snp', False):
            track_data = {
                "storeClass": "JBrowse/Store/SeqFeature/BAM",
                "urlTemplate": os.path.join('..', dest),
                "type": "JBrowse/View/Track/SNPCoverage",
                "key": key + " - SNPs/Coverage",
                "label": label + "_autosnp",
            }
            if 'category' in kwargs:
                track_data['category'] = kwargs['category']
            self._add_json(track_data)

    def add_vcf(self, data, key, **kwargs):
        label = hashlib.md5(data).hexdigest()
        dest = os.path.join('data', 'raw', os.path.basename(data))
        # ln?
        cmd = ['ln', '-s', data, dest]
        self.subprocess_check_call(cmd)
        cmd = ['bgzip', dest]
        self.subprocess_check_call(cmd)
        cmd = ['tabix', '-p', 'vcf', dest + '.gz']
        self.subprocess_check_call(cmd)

        track_data = {
            "key": key,
            "label": label,
            "urlTemplate": os.path.join('..', dest + '.gz'),
            "type": "JBrowse/View/Track/HTMLVariants",
            "storeClass": "JBrowse/Store/SeqFeature/VCFTabix",
        }
        track_data.update(kwargs)
        self._add_json(track_data)

    def add_features(self, data, key, format, **kwargs):
        label = hashlib.md5(data).hexdigest()
        cmd = ['perl', 'bin/flatfile-to-json.pl',
               TN_TABLE.get(format, 'gff'), data,
               '--trackLabel', label,
               '--key', key]

        config = {}
        if 'category' in kwargs:
            config['category'] = kwargs['category']

        if kwargs.get('match', False):
            clientConfig = {
                'label': 'description',
                'description': 'Hit_titles',
            }

            # Get min/max and build a scoring function since JBrowse doesn't
            min_val, max_val = self._min_max_gff(data)

            if min_val is not None and max_val is not None:
                MIN_MAX_OPACITY_MATH = Template("""
                var opacity = (score - ${min}) * (1/(${max} - ${min}));
                """).substitute({
                    'max': max_val,
                    'min': min_val,
                })

                red, green, blue = self._get_colours()
                color_function = COLOR_FUNCTION_TEMPLATE.substitute({
                    'score': "feature.get('score')",
                    'opacity': MIN_MAX_OPACITY_MATH,
                    'red': red,
                    'green': green,
                    'blue': blue,
                })

                clientConfig['color'] = color_function.replace('\n', '')

            config['glyph'] = 'JBrowse/View/FeatureGlyph/Segments'

            cmd += ['--clientConfig', json.dumps(clientConfig),
                    '--trackType', 'JBrowse/View/Track/CanvasFeatures'
                    ]

        cmd.extend(['--config', json.dumps(config)])

        self.subprocess_check_call(cmd)

    def process_annotations(self, data, key, format, **kwargs):
        if format in ('gff', 'gff3', 'bed'):
            self.add_features(data, key, format, **kwargs)
        elif format == 'bigwig':
            self.add_bigwig(data, key, **kwargs)
        elif format == 'bam':
            self.add_bam(data, key, **kwargs)
        elif format == 'blastxml':
            self.add_blastxml(data, key, **kwargs)
        elif format == 'vcf':
            self.add_vcf(data, key, **kwargs)

    def clone_jbrowse(self, jbrowse_dir, destination):
        # JBrowse seems to have included some bad symlinks, cp ignores bad symlinks
        # unlike copytree
        cmd = ['mkdir', '-p', destination]
        subprocess.check_call(cmd)
        cmd = ['cp', '-r', jbrowse_dir, destination]
        subprocess.check_call(cmd)
        cmd = ['mkdir', '-p', os.path.join(destination, 'JBrowse-1.11.6',
                                           'data', 'raw')]
        subprocess.check_call(cmd)

        # http://unix.stackexchange.com/a/38691/22785
        # JBrowse releases come with some broken symlinks
        cmd = ['find', destination, '-type', 'l', '-xtype', 'l', '-exec', 'rm', "'{}'", '+']
        subprocess.check_call(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument('genome', type=file, help='Input genome file')
    parser.add_argument('yaml', type=file, help='Track Configuration')

    parser.add_argument('--jbrowse', help='Folder containing a jbrowse release')
    parser.add_argument('--outdir', help='Output directory', default='out')
    args = parser.parse_args()

    jc = JbrowseConnector(
        jbrowse=args.jbrowse,
        jbrowse_dir=os.path.join(args.outdir, 'JBrowse-1.11.6'),
        outdir=args.outdir,
        genome=os.path.realpath(args.genome.name),
    )

    track_data = yaml.load(args.yaml)
    for track in track_data:
        path = os.path.realpath(track['file'])
        extra = track.get('options', {})
        if '__unused__' in extra:
            del extra['__unused__']

        for possible_partial_path in ('bam_index', 'parent'):
            if possible_partial_path in extra:
                extra[possible_partial_path] = os.path.realpath(
                    extra[possible_partial_path])
        extra['category'] = track.get('category', 'Default')

        jc.process_annotations(path, track['label'], track['ext'], **extra)

    print """
    <html>
        <body>
        <script type="text/javascript">
            window.location=JBrowse-1.11.6/index.html
        </script>
        <a href="JBrowse-1.11.6/index.html">Go to JBrowse</a>
        <p>Please note that JBrowse functions best on production Galaxy
        instances. The paste server used in development instances has issues
        serving the volumes of data regularly involved in JBrowse</p>
        </body>
    </html>
    """
