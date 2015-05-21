#!/usr/bin/env python
import os
import argparse
import subprocess
import hashlib
import tempfile
import json
import yaml

TN_TABLE = {
    'gff3': '--gff',
    'gff': '--gff',
    'bed': '--bed',
    # 'genbank': '--gbk',
}


def process_genome(jbrowse_dir, genome):
    subprocess.check_output(['perl', 'bin/prepare-refseqs.pl', '--fasta', genome], cwd=jbrowse_dir)


def _add_json(jbrowse_dir, json_data):
    if len(json_data.keys()) == 0:
        return

    tmp = tempfile.NamedTemporaryFile(delete=False)
    tmp.write(json.dumps(json_data))
    tmp.close()
    cmd = ['perl', 'bin/add-track-json.pl', tmp.name,
           os.path.join('data', 'trackList.json')]
    subprocess.check_call(cmd, cwd=jbrowse_dir)
    os.unlink(tmp.name)


def add_blastxml(jbrowse_dir, data, key, format, **kwargs):
    label = hashlib.md5(data).hexdigest()
    color_function = """
        function(feature, variableName, glyphObject, track) {
            var score = feature._parent.get('score');
            var opacity = 0;
            if(score == 0.0) {
                opacity = 1;
            } else{
                opacity = (20 - Math.log10(score)) / 180;
            }
            return 'rgba(109, 166, 166, ' + opacity + ')';
        }
    """
    clientConfig = {
        'label': 'description',
        'color': color_function.replace('\n', ''),
        'description': 'Hit_titles',
    }
    config = {
        'glyph': 'JBrowse/View/FeatureGlyph/Segments',
    }
    cmd = ['perl', 'bin/flatfile-to-json.pl',
           '--gff3', data,
           '--trackLabel', label,
           '--key', key,
           '--clientConfig', json.dumps(clientConfig),
           '--type', 'JBrowse/View/Track/CanvasFeatures',
           '--config', json.dumps(config),
           '--trackType', 'JBrowse/View/Track/CanvasFeatures'
           ]

    subprocess.check_call(cmd, cwd=jbrowse_dir)


def add_bigwig(jbrowse_dir, data, key, **kwargs):
    label = hashlib.md5(data).hexdigest()
    dest = os.path.join('data', 'raw', os.path.basename(data))
    # ln?
    cmd = ['cp', data, dest]
    subprocess.check_call(cmd, cwd=jbrowse_dir)

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

    _add_json(jbrowse_dir, track_data)


def add_bam(jbrowse_dir, data, key, **kwargs):
    label = hashlib.md5(data).hexdigest()
    dest = os.path.join('data', 'raw', os.path.basename(data))
    # ln?
    cmd = ['cp', data, dest]
    subprocess.check_call(cmd, cwd=jbrowse_dir)
    cmd = ['cp', data + '.bai', dest + '.bai']
    subprocess.check_call(cmd, cwd=jbrowse_dir)

    track_data = {
        "urlTemplate": os.path.join('..', dest),
        "key": key,
        "label": label,
        "type": "JBrowse/View/Track/Alignments2",
        "storeClass": "JBrowse/Store/SeqFeature/BAM",
    }
    _add_json(jbrowse_dir, track_data)

    if kwargs.get('auto_snp', False):
        track_data = {
            "storeClass": "JBrowse/Store/SeqFeature/BAM",
            "urlTemplate": os.path.join('..', dest),
            "type": "JBrowse/View/Track/SNPCoverage",
            "key": key + " - SNPs/Coverage",
            "label": label + "_autosnp",
        }
        _add_json(jbrowse_dir, track_data)


def add_vcf(jbrowse_dir, data, key, **kwargs):
    label = hashlib.md5(data).hexdigest()
    dest = os.path.join('data', 'raw', os.path.basename(data))
    # ln?
    cmd = ['cp', data, dest]
    subprocess.check_call(cmd, cwd=jbrowse_dir)
    cmd = ['bgzip', dest]
    subprocess.check_call(cmd, cwd=jbrowse_dir)
    cmd = ['tabix', '-p', 'vcf', dest + '.gz']
    subprocess.check_call(cmd, cwd=jbrowse_dir)

    track_data = {
        "key": key,
        "label": label,
        "urlTemplate": os.path.join('..', dest + '.gz'),
        "type": "JBrowse/View/Track/HTMLVariants",
        "storeClass": "JBrowse/Store/SeqFeature/VCFTabix",
    }
    track_data.update(kwargs)
    _add_json(jbrowse_dir, track_data)


def add_features(jbrowse_dir, data, key, format, **kwargs):
    label = hashlib.md5(data).hexdigest()
    cmd = ['perl', 'bin/flatfile-to-json.pl', TN_TABLE.get(format, 'gff'),
           data, '--trackLabel', label, '--key', key]
    subprocess.check_call(cmd, cwd=jbrowse_dir)


def process_annotations(jbrowse_dir, data, key, format,
                        **kwargs):
    if format in ('gff', 'gff3', 'bed'):
        add_features(jbrowse_dir, data, key, format, **kwargs)
    elif format == 'bigwig':
        add_bigwig(jbrowse_dir, data, key, **kwargs)
    elif format == 'bam':
        add_bam(jbrowse_dir, data, key, **kwargs)
    elif format in ('blastxml', ):
        add_blastxml(jbrowse_dir, data, key, **kwargs)
    elif format == 'vcf':
        add_vcf(jbrowse_dir, data, key, **kwargs)
        pass


def clone_jbrowse(jbrowse_dir, destination):
    # JBrowse seems to have included some bad symlinks, cp ignores bad symlinks
    # unlike copytree
    cmd = ['mkdir', '-p', destination]
    subprocess.check_output(cmd)
    cmd = ['cp', '-r', jbrowse_dir, destination]
    subprocess.check_output(cmd)
    cmd = ['mkdir', '-p', os.path.join(destination, 'JBrowse-1.11.6', 'data',
                                       'raw')]
    subprocess.check_output(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument('genome', type=file, help='Input genome file')
    parser.add_argument('yaml', type=file, help='Track Configuration')

    parser.add_argument('--jbrowse', help='Folder containing a jbrowse release')
    parser.add_argument('--outdir', help='Output directory', default='out')
    args = parser.parse_args()

    jbrowse_dir = os.path.join(args.outdir, 'JBrowse-1.11.6')
    clone_jbrowse(args.jbrowse, args.outdir)

    process_genome(jbrowse_dir, os.path.realpath(args.genome.name))

    track_data = yaml.load(args.yaml)
    for track in track_data:
        path = os.path.realpath(track['file'])
        extra = track.get('options', {})
        process_annotations(jbrowse_dir, path, track['label'], track['ext'], **extra)

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
