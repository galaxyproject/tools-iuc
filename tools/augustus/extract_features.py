#!/usr/bin/env python

import argparse
import sys
import textwrap


def main( args ):
    """
    Extract the protein and coding section from an augustus gff, gtf file
    Example file:
HS04636	AUGUSTUS	stop_codon	6901	6903	.	+	0	Parent=g1.t1
HS04636	AUGUSTUS	transcription_end_site	8857	8857	.	+	.	Parent=g1.t1
# protein sequence = [MLARALLLCAVLALSHTANPCCSHPCQNRGVCMSVGFDQYKCDCTRTGFYGENCSTPEFLTRIKLFLKPTPNTVHYIL
# THFKGFWNVVNNIPFLRNAIMSYVLTSRSHLIDSPPTYNADYGYKSWEAFSNLSYYTRALPPVPDDCPTPLGVKGKKQLPDSNEIVEKLLLRRKFIPD
# PQGSNMMFAFFAQHFTHQFFKTDHKRGPAFTNGLGHGVDLNHIYGETLARQRKLRLFKDGKMKYQIIDGEMYPPTVKDTQAEMIYPPQVPEHLRFAVG
# QEVFGLVPGLMMYATIWLREHNRVCDVLKQEHPEWGDEQLFQTSRLILIGETIKIVIEDYVQHLSGYHFKLKFDPELLFNKQFQYQNRIAAEFNTLYH
# WHPLLPDTFQIHDQKYNYQQFIYNNSILLEHGITQFVESFTRQIAGRVAGGRNVPPAVQKVSQASIDQSRQMKYQSFNEYRKRFMLKPYESFEELTGE
# KEMSAELEALYGDIDAVELYPALLVEKPRPDAIFGETMVEVGAPFSLKGLMGNVICSPAYWKPSTFGGEVGFQIINTASIQSLICNNVKGCPFTSFSV
# PDPELIKTVTINASSSRSGLDDINPTVLLKERSTEL]
# end gene g1
###
#
# ----- prediction on sequence number 2 (length = 2344, name = HS08198) -----
#
# Predicted genes for sequence number 2 on both strands
# start gene g2
HS08198	AUGUSTUS	gene	86	2344	1	+	.	ID=g2
HS08198	AUGUSTUS	transcript	86	2344	.	+	.	ID=g2.t1;Parent=g2
HS08198	AUGUSTUS	transcription_start_site	86	86	.	+	.	Parent=g2.t1
HS08198	AUGUSTUS	exon	86	582	.	+	.	Parent=g2.t1
HS08198	AUGUSTUS	start_codon	445	447	.	+	0	Parent=g2.t1
    """
    protein_seq = ''
    coding_seq = ''
    if args.protein:
        po = open( args.protein, 'w+' )
    if args.codingseq:
        co = open( args.codingseq, 'w+' )

    for line in sys.stdin:
        # protein- and coding-sequence are stored as comments
        if line.startswith('#'):
            line = line[2:].strip()
            if line.startswith('start gene'):
                gene_name = line[11:].strip()

            if protein_seq:
                if line.endswith(']'):
                    protein_seq += line[:-1]
                    po.write( '>%s\n%s\n' % (gene_name, '\n'.join( textwrap.wrap( protein_seq, 80 ) ) ) )
                    protein_seq = ''
                else:
                    protein_seq += line

            if coding_seq:
                if line.endswith(']'):
                    coding_seq += line[:-1]
                    co.write( '>%s\n%s\n' % (gene_name, '\n'.join( textwrap.wrap( coding_seq, 80 ) ) ) )
                    coding_seq = ''
                else:
                    coding_seq += line

            if args.protein and line.startswith('protein sequence = ['):
                if line.endswith(']'):
                    protein_seq = line[20:-1]
                    po.write( '>%s\n%s\n' % (gene_name, '\n'.join( textwrap.wrap( protein_seq, 80 ) ) ) )
                    protein_seq = ''
                else:
                    line = line[20:]
                    protein_seq = line

            if args.codingseq and line.startswith('coding sequence = ['):
                if line.endswith(']'):
                    coding_seq = line[19:-1]
                    co.write( '>%s\n%s\n' % (gene_name, '\n'.join( textwrap.wrap( coding_seq, 80 ) ) ) )
                    coding_seq = ''
                else:
                    line = line[19:]
                    coding_seq = line

    if args.codingseq:
        co.close()
    if args.protein:
        po.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--protein', help='Path to the protein file.')
    parser.add_argument('-c', '--codingseq', help='Path to the coding file.')

    args = parser.parse_args()
    main( args )
