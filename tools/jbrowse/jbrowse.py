#!/usr/bin/env python
import os
import argparse
import subprocess
import hashlib
import tempfile
import json

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

    print json_data
    tmp = tempfile.NamedTemporaryFile(delete=False)
    tmp.write(json.dumps(json_data))
    tmp.close()
    cmd = ['perl', 'bin/add-track-json.pl', tmp.name,
           os.path.join('data', 'trackList.json')]
    subprocess.check_call(cmd, cwd=jbrowse_dir)
    os.unlink(tmp.name)


def add_bigwig(jbrowse_dir, data, key, **kwargs):
    label = hashlib.md5(data).hexdigest()
    source = data
    dest = os.path.join('data', 'raw', os.path.basename(data))
    # ln?
    cmd = ['cp', source, dest]
    subprocess.check_call(cmd, cwd=jbrowse_dir)

    track_data = {
        "autoscale": "local",
        "label": label,
        "urlTemplate": os.path.join('..', dest),
        "bicolor_pivot": "zero",
        "storeClass": "JBrowse/Store/SeqFeature/BigWig",
        "type": "JBrowse/View/Track/Wiggle/Density",
        "key": "bigwig"
    }
    track_data.update(kwargs)
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
    elif format in ('bigwig', 'wig'):
        add_bigwig(jbrowse_dir, data, key, **kwargs)
    elif format in ('bam', ):
        add_bam(jbrowse_dir, data, key, **kwargs)
    elif format in ('blastxml', ):
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

    parser.add_argument('--annotation', type=file, nargs='*',
                        help='annotation data')
    parser.add_argument('--annotation_format', choices=['gff3', 'gff', 'bed', 'bam', 'wig', 'bigwig'],
                        nargs='*', help='annotation format')
    parser.add_argument('--annotation_label', type=str, nargs='*',
                        help='annotation label')
    parser.add_argument('--annotation_extra', type=str, nargs='*',
                        help='annotation extra data')

    parser.add_argument('--jbrowse', help='Folder containing a jbrowse release')
    parser.add_argument('--outdir', help='Output directory', default='out')
    args = parser.parse_args()

    jbrowse_dir = os.path.join(args.outdir, 'JBrowse-1.11.6')
    clone_jbrowse(args.jbrowse, args.outdir)

    process_genome(jbrowse_dir, os.path.realpath(args.genome.name))

    for annotation, format, label, extra in zip(args.annotation,
                                                args.annotation_format,
                                                args.annotation_label,
                                                args.annotation_extra):
        path = os.path.realpath(annotation.name)

        # Take in json formatted "extra" that goes straight into trackData.json
        real_extra = {}
        if extra != "None":
            real_extra = json.loads(extra)

        process_annotations(jbrowse_dir, path, label, format, **real_extra)

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
