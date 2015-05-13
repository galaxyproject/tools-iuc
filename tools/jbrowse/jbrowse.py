#!/usr/bin/env python
import os
import argparse
import subprocess
import hashlib

TN_TABLE = {
    'gff3': '--gff',
    'gff': '--gff',
    'bed': '--bed',
    # 'genbank': '--gbk',
}


def process_genome(jbrowse_dir, genome):
    subprocess.check_output(['perl', 'bin/prepare-refseqs.pl', '--fasta', genome], cwd=jbrowse_dir)


def process_annotations(jbrowse_dir, data, label, format,
                        **kwargs):
    key = hashlib.md5(data).hexdigest()

    if format in ('gff', 'gff3', 'bed'):
        cmd = ['perl', 'bin/flatfile-to-json.pl', TN_TABLE.get(format, 'gff'),
               data, '--trackLabel', key, '--key', label]
    elif format in ('bigwig', 'wig'):
        cmd = []
    elif format in ('bam', ):
        cmd = []
    elif format in ('blastxml', ):
        cmd = []
    subprocess.check_output(cmd, cwd=jbrowse_dir)


def clone_jbrowse(jbrowse_dir, destination):
    # JBrowse seems to have included some bad symlinks, cp ignores bad symlinks
    # unlike copytree
    cmd = ['mkdir', '-p', destination]
    subprocess.check_output(cmd)
    cmd = ['cp', '-r', jbrowse_dir, destination]
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

    parser.add_argument('--jbrowse', help='Folder containing a jbrowse release')
    parser.add_argument('--outdir', help='Output directory', default='out')
    args = parser.parse_args()

    jbrowse_dir = os.path.join(args.outdir, 'JBrowse-1.11.6')
    clone_jbrowse(args.jbrowse, args.outdir)

    process_genome(jbrowse_dir, os.path.realpath(args.genome.name))

    for annotation, format, label in zip(args.annotation, args.annotation_format, args.annotation_label):
        path = os.path.realpath(annotation.name)
        process_annotations(jbrowse_dir, path, label, format)

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
