#!/usr/bin/env python
import os
import argparse
import subprocess
import hashlib
import tempfile
import json
import yaml

COLOR_FUNCTION_TEMPLATE = """
function(feature, variableName, glyphObject, track) {
    var score = %s;
    %s
    return 'rgba(%s, %s, %s, ' + opacity + ')';
}
"""

BLAST_OPACITY_MATH = """
var opacity = 0;
if(score == 0.0) {
    opacity = 1;
} else{
    opacity = (20 - Math.log10(score)) / 180;
}
"""

CANVAS_ARGS = [
    '--config', json.dumps({'glyph': 'JBrowse/View/FeatureGlyph/Segments'}),
    '--type', 'JBrowse/View/Track/CanvasFeatures',
    '--trackType', 'JBrowse/View/Track/CanvasFeatures'
]

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

    def _get_colours(self):
        r, g, b = BREWER_COLOUR_SCHEMES[self.brewer_colour_idx]
        self.brewer_colour_idx += 1
        return r, g, b

    def process_genome(self):
        subprocess.check_output(['perl', 'bin/prepare-refseqs.pl', '--fasta',
                                 self.genome_path], cwd=self.jbrowse_dir)

    def _add_json(self, json_data):
        if len(json_data.keys()) == 0:
            return

        tmp = tempfile.NamedTemporaryFile(delete=False)
        tmp.write(json.dumps(json_data))
        tmp.close()
        cmd = ['perl', 'bin/add-track-json.pl', tmp.name,
               os.path.join('data', 'trackList.json')]
        subprocess.check_call(cmd, cwd=self.jbrowse_dir)
        os.unlink(tmp.name)

    def add_blastxml(self, data, key, format, **kwargs):
        gff3_unrebased = tempfile.NamedTemporaryFile(delete=False)
        cmd = ['python', os.path.join(INSTALLED_TO, 'blastxml_to_gapped_gff3.py'),
               '--trim_end', '--min_gap', '10', data]
        subprocess.check_call(cmd, cwd=self.jbrowse_dir, stdout=gff3_unrebased)
        gff3_unrebased.close()

        gff3_rebased = tempfile.NamedTemporaryFile(delete=False)
        cmd = ['python', os.path.join(INSTALLED_TO, 'gff3_rebase.py')]
        if kwargs['protein']:
            cmd.append('--protein2dna')
        cmd.extend(['--trim_end', '--min_gap', kwargs['min_gap'], kwargs['parent'], gff3_unrebased.name])
        subprocess.check_call(cmd, cwd=self.jbrowse_dir, stdout=gff3_rebased)
        gff3_rebased.close()

        label = hashlib.md5(data).hexdigest()

        color_function = COLOR_FUNCTION_TEMPLATE % \
            ["feature._parent.get('score')",
             BLAST_OPACITY_MATH].extend(self._get_colours())

        clientConfig = {
            'label': 'description',
            'color': color_function.replace('\n', ''),
            'description': 'Hit_titles',
        }
        cmd = ['perl', 'bin/flatfile-to-json.pl',
               '--gff3', gff3_rebased.name,
               '--trackLabel', label,
               '--key', key,
               '--clientConfig', json.dumps(clientConfig),
               ] + CANVAS_ARGS

        subprocess.check_call(cmd, cwd=self.jbrowse_dir)

    def add_bigwig(self, data, key, **kwargs):
        label = hashlib.md5(data).hexdigest()
        dest = os.path.join('data', 'raw', os.path.basename(data))
        # ln?
        cmd = ['cp', data, dest]
        subprocess.check_call(cmd, cwd=self.jbrowse_dir)

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
        cmd = ['cp', data, dest]
        subprocess.check_call(cmd, cwd=self.jbrowse_dir)

        bai_source = kwargs['bam_index']
        cmd = ['cp', bai_source, dest + '.bai']
        subprocess.check_call(cmd, cwd=self.jbrowse_dir)

        track_data = {
            "urlTemplate": os.path.join('..', dest),
            "key": key,
            "label": label,
            "type": "JBrowse/View/Track/Alignments2",
            "storeClass": "JBrowse/Store/SeqFeature/BAM",
        }
        self._add_json(track_data)

        if kwargs.get('auto_snp', False):
            track_data = {
                "storeClass": "JBrowse/Store/SeqFeature/BAM",
                "urlTemplate": os.path.join('..', dest),
                "type": "JBrowse/View/Track/SNPCoverage",
                "key": key + " - SNPs/Coverage",
                "label": label + "_autosnp",
            }
            self._add_json(track_data)

    def add_vcf(self, data, key, **kwargs):
        label = hashlib.md5(data).hexdigest()
        dest = os.path.join('data', 'raw', os.path.basename(data))
        # ln?
        cmd = ['cp', data, dest]
        subprocess.check_call(cmd, cwd=self.jbrowse_dir)
        cmd = ['bgzip', dest]
        subprocess.check_call(cmd, cwd=self.jbrowse_dir)
        cmd = ['tabix', '-p', 'vcf', dest + '.gz']
        subprocess.check_call(cmd, cwd=self.jbrowse_dir)

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
        cmd = ['perl', 'bin/flatfile-to-json.pl', TN_TABLE.get(format, 'gff'),
               data, '--trackLabel', label, '--key', key]
        subprocess.check_call(cmd, cwd=self.jbrowse_dir)

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
        subprocess.check_output(cmd)
        cmd = ['cp', '-r', jbrowse_dir, destination]
        subprocess.check_output(cmd)
        cmd = ['mkdir', '-p', os.path.join(destination, 'JBrowse-1.11.6',
                                           'data', 'raw')]
        subprocess.check_output(cmd)

        # http://unix.stackexchange.com/a/38691/22785
        # JBrowse releases come with some broken symlinks
        cmd = ['find', destination, '-type', 'l', '-xtype', 'l', '-exec', 'rm', "'{}'", '\;']
        subprocess.check_output(cmd)


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
        jc.process_annotations(path, track['label'], track['ext'], **extra)

    print """
    <html>
        <body>
        <script type="text/javascript">
            window.location=JBrowse-1.11.6/index.html
        </script>
        <a href="JBrowse-1.11.6/index.html">Go to JBrowse</a>
        </body>
    </html>
    """
