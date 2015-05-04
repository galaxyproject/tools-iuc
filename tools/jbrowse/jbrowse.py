#!/usr/bin/env python
import os
import shutil
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


def process_annotations(jbrowse_dir, annot_file, annot_label, annot_format,
                        **kwargs):
    key = hashlib.md5(annot_file).hexdigest()

    cmd = ['perl', 'bin/flatfile-to-json.pl', TN_TABLE.get(annot_format, 'gff'),
           annot_file, '--trackLabel', key, '--key', annot_label]
    subprocess.check_output(cmd, cwd=jbrowse_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument('genome', type=file, help='Input genome file')

    parser.add_argument('--gff3', type=file, nargs='*', help='GFF3/BED/GBK File')
    parser.add_argument('--gff3_format', choices=['gff3', 'gff', 'bed', 'gbk'], nargs='*', help='GFF3/BED/GBK Format')
    parser.add_argument('--gff3_label', type=str, nargs='*', help='GFF3/BED/GBK Label')

    parser.add_argument('--jbrowse', help='Folder containing a jbrowse release')
    parser.add_argument('--outdir', help='Output directory', default='out')
    args = parser.parse_args()

    jbrowse_dir = os.path.join(args.outdir, 'JBrowse-1.11.6')
    shutil.copytree(args.jbrowse, jbrowse_dir)

    process_genome(jbrowse_dir, os.path.realpath(args.genome.name))

    for gff3, gff3_format, gff3_label in zip(args.gff3, args.gff3_format, args.gff3_label):
        gff3_path = os.path.realpath(gff3.name)
        process_annotations(jbrowse_dir, gff3_path, gff3_label, gff3_format)

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
