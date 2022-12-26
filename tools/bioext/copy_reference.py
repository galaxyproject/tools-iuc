import argparse
import os
import shutil

from BioExt.references import cov2, hxb2, nl4_3


references = {
    'HXB2_env': hxb2.env,
    'HXB2_gag': hxb2.gag,
    'HXB2_int': hxb2.int,
    'HXB2_nef': hxb2.nef,
    'HXB2_pol': hxb2.pol,
    'HXB2_pr': hxb2.pr,
    'HXB2_prrt': hxb2.prrt,
    'HXB2_rev': hxb2.rev,
    'HXB2_rt': hxb2.rt,
    'HXB2_tat': hxb2.tat,
    'HXB2_vif': hxb2.vif,
    'HXB2_vpr': hxb2.vpr,
    'HXB2_vpu': hxb2.vpu,
    'NL4-3_prrt': nl4_3.prrt,
    'CoV2-3C': cov2.threeC,
    'CoV2-E': cov2.E,
    'CoV2-endornase': cov2.endornase,
    'CoV2-exonuclease': cov2.exonuclease,
    'CoV2-helicase': cov2.helicase,
    'CoV2-leader': cov2.leader,
    'CoV2-methyltransferase': cov2.methyltransferase,
    'CoV2-M': cov2.M,
    'CoV2-N': cov2.N,
    'CoV2-nsp10': cov2.nsp10,
    'CoV2-nsp2': cov2.nsp2,
    'CoV2-nsp3': cov2.nsp3,
    'CoV2-nsp4': cov2.nsp4,
    'CoV2-nsp6': cov2.nsp6,
    'CoV2-nsp7': cov2.nsp7,
    'CoV2-nsp8': cov2.nsp8,
    'CoV2-nsp9': cov2.nsp9,
    'CoV2-ORF10': cov2.ORF10,
    'CoV2-ORF1a': cov2.ORF1a,
    'CoV2-ORF1b': cov2.ORF1b,
    'CoV2-ORF3a': cov2.ORF3a,
    'CoV2-ORF5': cov2.ORF5,
    'CoV2-ORF6': cov2.ORF6,
    'CoV2-ORF7a': cov2.ORF7a,
    'CoV2-ORF7b': cov2.ORF7b,
    'CoV2-ORF8': cov2.ORF8,
    'CoV2-RdRp': cov2.RdRp,
    'CoV2-S': cov2.S
}

parser = argparse.ArgumentParser()
parser.add_argument('--reference', dest='reference', action='store', type=str)
parser.add_argument('--dataset', dest='dataset', action='store', type=str)
args = parser.parse_args()
reference = os.path.abspath(references[args.reference]._seqpath)

shutil.copy(reference, os.path.abspath(args.dataset))
