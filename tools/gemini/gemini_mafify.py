import string
import sys


so_to_maf = {
    'splice_acceptor_variant': 'Splice_Site',
    'splice_donor_variant': 'Splice_Site',
    'transcript_ablation': 'Splice_Site',
    'exon_loss_variant': 'Splice_Site',
    'stop_gained': 'Nonsense_Mutation',
    'stop_lost': 'Nonstop_Mutation',
    'frameshift_variant': 'Frame_Shift_',
    'initiator_codon_variant': 'Translation_Start_Site',
    'start_lost': 'Translation_Start_Site',
    'inframe_insertion': 'In_Frame_Ins',
    'inframe_deletion': 'In_Frame_Del',
    'conservative_inframe_insertion': 'In_Frame_Ins',
    'conservative_inframe_deletion': 'In_Frame_Del',
    'disruptive_inframe_insertion': 'In_Frame_Ins',
    'disruptive_inframe_deletion': 'In_Frame_Del',
    'missense_variant': 'Missense_Mutation',
    'coding_sequence_variant': 'Missense_Mutation',
    'conservative_missense_variant': 'Missense_Mutation',
    'rare_amino_acid_variant': 'Missense_Mutation',
    'transcript_amplification': 'Intron',
    'intron_variant': 'Intron',
    'INTRAGENIC': 'Intron',
    'intragenic_variant': 'Intron',
    'splice_region_variant': 'Splice_Region',
    'mature_miRNA_variant': 'RNA',
    'exon_variant': 'RNA',
    'non_coding_exon_variant': 'RNA',
    'non_coding_transcript_exon_variant': 'RNA',
    'non_coding_transcript_variant': 'RNA',
    'nc_transcript_variant': 'RNA',
    'stop_retained_variant': 'Silent',
    'synonymous_variant': 'Silent',
    'NMD_transcript_variant': 'Silent',
    'incomplete_terminal_codon_variant': 'Silent',
    '5_prime_UTR_variant': "5'UTR",
    '5_prime_UTR_premature_start_codon_gain_variant': "5'UTR",
    '3_prime_UTR_variant': "3'UTR",
    'intergenic_variant': 'IGR',
    'intergenic_region': 'IGR',
    'regulatory_region_variant': 'IGR',
    'regulatory_region': 'IGR',
    'TF_binding_site_variant': 'IGR',
    'upstream_gene_variant': "5'Flank",
    'downstream_gene_variant': "3'Flank",
}


class VariantEffect():
    def __init__(self, variant_type):
        self.variant_type = variant_type.capitalize()
        assert self.variant_type in ['Snp', 'Ins', 'Del']

    def __getitem__(self, so_effect):
        if so_effect not in so_to_maf or (
            'frame' in so_effect and self.variant_type == 'Snp'
        ):
            return 'Targeted_Region'

        ret = so_to_maf[so_effect]
        if ret == 'Frame_Shift_':
            ret += self.variant_type
        return ret


infile = sys.argv[1]
if len(sys.argv) > 2:
    tumor_sample_name = sys.argv[2]
if len(sys.argv) > 3:
    normal_sample_name = sys.argv[3]

start_pos_idx = None
ref_idx = None
alt_idx = None
variant_type_idx = None
variant_classification_idx = None
gt_alt_depths_idx = {}
gt_ref_depths_idx = {}
gts_idx = {}
samples = set()
required_fields = [
    'Hugo_Symbol',
    'NCBI_Build',
    'Variant_Type',
    'Variant_Classification',
    'Tumor_Sample_Barcode',
    'HGVSp_Short'
]


with open(infile) as data_in:
    cols = data_in.readline().rstrip().split('\t')
    for field in required_fields:
        if field not in cols:
            raise IndexError(
                'Cannot generate valid MAF without the following input '
                'columns: {0}.\n'
                'Missing column: "{1}"'
                .format(required_fields, field)
            )
    for i, col in enumerate(cols):
        if col == 'Variant_Type':
            variant_type_idx = i
        elif col == 'Variant_Classification':
            variant_classification_idx = i
        elif col == 'Start_Position':
            start_pos_idx = i
        elif col == 'Reference_Allele':
            ref_idx = i
        elif col == 'alt':
            alt_idx = i
        else:
            column, _, sample = col.partition('.')
            if sample:
                if column == 'gt_alt_depths':
                    gt_alt_depths_idx[sample] = i
                elif column == 'gt_ref_depths':
                    gt_ref_depths_idx[sample] = i
                elif column == 'gts':
                    gts_idx[sample] = i
                else:
                    # not a recognized sample-specific column
                    continue
                samples.add(sample)

    if ref_idx is None:
        raise IndexError('Input file does not have a column "Reference_Allele".')

    if not tumor_sample_name:
        if normal_sample_name:
            raise ValueError(
                'Normal sample name requires the tumor sample name to be '
                'specified, too.'
            )
        if len(samples) > 1:
            raise ValueError(
                'A tumor sample name is required with more than one sample '
                'in the input.'
            )
        if samples:
            # There is a single sample with genotype data.
            # Assume its the tumor sample.
            tumor_sample_name = next(iter(samples))
    else:
        if tumor_sample_name not in samples:
            raise ValueError(
                'Could not find information about the specified tumor sample '
                'in the input.'
            )
        if tumor_sample_name == normal_sample_name:
            raise ValueError(
                'Need different names for the normal and the tumor sample.'
            )

    if normal_sample_name and normal_sample_name not in samples:
        raise ValueError(
            'Could not find information about the specified normal sample '
            'in the input.'
        )

    # All input data checks passed!
    # Now extract just the relevant index numbers for the tumor/normal pair
    gts_idx = (
        gts_idx.get(tumor_sample_name, alt_idx),
        gts_idx.get(normal_sample_name)
    )
    gt_alt_depths_idx = (
        gt_alt_depths_idx.get(tumor_sample_name),
        gt_alt_depths_idx.get(normal_sample_name)
    )
    gt_ref_depths_idx = (
        gt_ref_depths_idx.get(tumor_sample_name),
        gt_ref_depths_idx.get(normal_sample_name)
    )

    # Echo all MAF column names
    cols_to_print = []
    for n in range(len(cols)):
        if n in gts_idx:
            continue
        if n in gt_alt_depths_idx:
            continue
        if n in gt_ref_depths_idx:
            continue
        if n != alt_idx:
            cols_to_print.append(n)

    print('\t'.join([cols[n] for n in cols_to_print]))

    for line in data_in:
        cols = line.rstrip().split('\t')

        gt_alt_depths = [
            int(cols[ad_idx]) if ad_idx else ''
            for ad_idx in gt_alt_depths_idx
        ]
        gt_ref_depths = [
            int(cols[rd_idx]) if rd_idx else ''
            for rd_idx in gt_ref_depths_idx
        ]

        gts = [
            ['', ''],
            ['', '']
        ]
        for n, gt_idx in enumerate(gts_idx):
            if gt_idx:
                gt_sep = '/' if '/' in cols[gt_idx] else '|'
                allele1, _, allele2 = [
                    '' if allele == '.' else allele
                    for allele in cols[gt_idx].partition(gt_sep)
                ]
                # follow cBioportal recommendation to leave allele1 empty
                # when information is not avaliable
                if not allele2:
                    gts[n] = [allele2, allele1]
                else:
                    gts[n] = [allele1, allele2]
        if not gts:
            gts = [['', ''], ['', '']]

        if cols[variant_type_idx].lower() in ['ins', 'del']:
            # transform VCF-style indel representations into MAF ones
            ref_allele = cols[ref_idx]
            for n, nucs in enumerate(
                zip(
                    ref_allele,
                    *[allele for gt in gts for allele in gt if allele]
                )
            ):
                if any(nuc != nucs[0] for nuc in nucs[1:]):
                    break
            else:
                n += 1
            if n > 0:
                cols[ref_idx] = cols[ref_idx][n:] or '-'
                for gt in gts:
                    for idx, allele in enumerate(gt):
                        if allele:
                            gt[idx] = allele[n:] or '-'
                if cols[ref_idx] == '-':
                    n -= 1
                cols[start_pos_idx] = str(int(cols[start_pos_idx]) + n)

        # in-place substitution of so_effect with MAF effect
        cols[variant_classification_idx] = VariantEffect(
            cols[variant_type_idx]
        )[cols[variant_classification_idx]]
        ret_line = '\t'.join([cols[n] for n in cols_to_print])

        field_formatters = {
            'tumor_seq_allele1': gts[0][0],
            'tumor_seq_allele2': gts[0][1],
            'match_norm_seq_allele1': gts[1][0],
            'match_norm_seq_allele2': gts[1][1],
            't_alt_count': gt_alt_depths[0],
            'n_alt_count': gt_alt_depths[1],
            't_ref_count': gt_ref_depths[0],
            'n_ref_count': gt_ref_depths[1],
        }

        print(
            # use safe_substitute here to avoid key errors with column content
            # looking like unknown placeholders
            string.Template(ret_line).safe_substitute(field_formatters)
        )
