#!/usr/bin/env python

import optparse
import os
import shutil
import subprocess
import sys
import tempfile


def parse_gff_attributes( attr_str ):
    """
    Parses a GFF/GTF attribute string and returns a dictionary of name-value
    pairs. The general format for a GFF3 attributes string is

        name1=value1;name2=value2

    The general format for a GTF attribute string is

        name1 "value1" ; name2 "value2"

    The general format for a GFF attribute string is a single string that
    denotes the interval's group; in this case, method returns a dictionary
    with a single key-value pair, and key name is 'group'
    """
    attributes_list = attr_str.split(";")
    attributes = {}
    for name_value_pair in attributes_list:
        # Try splitting by '=' (GFF3) first because spaces are allowed in GFF3
        # attribute; next, try double quotes for GTF.
        pair = name_value_pair.strip().split("=")
        if len( pair ) == 1:
            pair = name_value_pair.strip().split("\"")
        if len( pair ) == 1:
            # Could not split for some reason -- raise exception?
            continue
        if pair == '':
            continue
        name = pair[0].strip()
        if name == '':
            continue
        # Need to strip double quote from values
        value = pair[1].strip(" \"")
        attributes[ name ] = value

    if len( attributes ) == 0:
        # Could not split attributes string, so entire string must be
        # 'group' attribute. This is the case for strictly GFF files.
        attributes['group'] = attr_str
    return attributes


def gff_attributes_to_str( attrs, gff_format ):
    """
    Convert GFF attributes to string. Supported formats are GFF3, GTF.
    """
    if gff_format == 'GTF':
        format_string = '%s "%s"'
        # Convert group (GFF) and ID, parent (GFF3) attributes to transcript_id, gene_id
        id_attr = None
        if 'group' in attrs:
            id_attr = 'group'
        elif 'ID' in attrs:
            id_attr = 'ID'
        elif 'Parent' in attrs:
            id_attr = 'Parent'
        if id_attr:
            attrs['transcript_id'] = attrs['gene_id'] = attrs[id_attr]
    elif gff_format == 'GFF3':
        format_string = '%s=%s'
    attrs_strs = []
    for name, value in attrs.items():
        attrs_strs.append( format_string % ( name, value ) )
    return " ; ".join( attrs_strs )


def stop_err(msg):
    sys.exit("%s\n" % msg)


def __main__():
    parser = optparse.OptionParser()
    parser.add_option('-1', '--input', dest='input', help=' file of RNA-Seq read alignments in the SAM format. SAM is a standard short read alignment, that allows aligners to attach custom tags to individual alignments, and Cufflinks requires that the alignments you supply have some of these tags. Please see Input formats for more details.')
    parser.add_option('-I', '--max-intron-length', dest='max_intron_len', help='The minimum intron length. Cufflinks will not report transcripts with introns longer than this, and will ignore SAM alignments with REF_SKIP CIGAR operations longer than this. The default is 300,000.')
    parser.add_option('-F', '--min-isoform-fraction', dest='min_isoform_fraction', help='After calculating isoform abundance for a gene, Cufflinks filters out transcripts that it believes are very low abundance, because isoforms expressed at extremely low levels often cannot reliably be assembled, and may even be artifacts of incompletely spliced precursors of processed transcripts. This parameter is also used to filter out introns that have far fewer spliced alignments supporting them. The default is 0.05, or 5% of the most abundant isoform (the major isoform) of the gene.')
    parser.add_option('-j', '--pre-mrna-fraction', dest='pre_mrna_fraction', help='Some RNA-Seq protocols produce a significant amount of reads that originate from incompletely spliced transcripts, and these reads can confound the assembly of fully spliced mRNAs. Cufflinks uses this parameter to filter out alignments that lie within the intronic intervals implied by the spliced alignments. The minimum depth of coverage in the intronic region covered by the alignment is divided by the number of spliced reads, and if the result is lower than this parameter value, the intronic alignments are ignored. The default is 5%.')
    parser.add_option('-p', '--num-threads', dest='num_threads', help='Use this many threads to align reads. The default is 1.')
    parser.add_option('-G', '--GTF', dest='GTF', help='Tells Cufflinks to use the supplied reference annotation to estimate isoform expression. It will not assemble novel transcripts, and the program will ignore alignments not structurally compatible with any reference transcript.')
    parser.add_option("--compatible-hits-norm", dest='compatible_hits_norm', action="store_true", help='Count hits compatible with reference RNAs only')
    parser.add_option('-g', '--GTF-guide', dest='GTFguide', help='use reference transcript annotation to guide assembly')
    parser.add_option("--3-overhang-tolerance", dest='three_overhang_tolerance', help='The number of bp allowed to overhang the 3prime end of a reference transcript when determining if an assembled transcript should be merged with it (ie, the assembled transcript is not novel). The default is 600 bp.')
    parser.add_option("--intron-overhang-tolerance", dest='intron_overhang_tolerance', help='The number of bp allowed to enter the intron of a reference transcript when determining if an assembled transcript should be merged with it (ie, the assembled transcript is not novel). The default is 50 bp.')
    parser.add_option("--no-faux-reads", dest='no_faux_reads', help='This option disables tiling of the reference transcripts with faux reads. Use this if you only want to use sequencing reads in assembly but do not want to output assembled transcripts that lay within reference transcripts. All reference transcripts in the input annotation will also be included in the output.')
    parser.add_option('-u', '--multi-read-correct', dest='multi_read_correct', action="store_true", help='Tells Cufflinks to do an initial estimation procedure to more accurately weight reads mapping to multiple locations in the genome')

    # Normalization options.
    parser.add_option("--no-effective-length-correction", dest="no_effective_length_correction", action="store_true")
    parser.add_option("--no-length-correction", dest="no_length_correction", action="store_true")

    # Wrapper / Galaxy options.
    parser.add_option('-A', '--assembled-isoforms-output', dest='assembled_isoforms_output_file', help='Assembled isoforms output file; formate is GTF.')

    # Advanced Options:
    parser.add_option("--library-type", dest="library_type", help=' library prep used for input reads, default fr-unstranded')
    parser.add_option('-M', '--mask-file', dest='mask_file', help='Tells Cufflinks to ignore all reads that could have come from transcripts in this GTF file. \
                                                                   We recommend including any annotated rRNA, mitochondrial transcripts other abundant transcripts \
                                                                   you wish to ignore in your analysis in this file. Due to variable efficiency of mRNA enrichment \
                                                                   methods and rRNA depletion kits, masking these transcripts often improves the overall robustness \
                                                                   of transcript abundance estimates.')
    parser.add_option('-m', '--inner-mean-dist', dest='inner_mean_dist', help='This is the expected (mean) inner distance between mate pairs. \
                                                                                For, example, for paired end runs with fragments selected at 300bp, \
                                                                                where each end is 50bp, you should set -r to be 200. The default is 45bp.')  # cufflinks: --frag-len-mean

    parser.add_option('-s', '--inner-dist-std-dev', dest='inner_dist_std_dev', help='The standard deviation for the distribution on inner distances between mate pairs. The default is 20bp.')  # cufflinks: --frag-len-std-dev
    parser.add_option('--max-mle-iterations', dest='max_mle_iterations', help='Sets the number of iterations allowed during maximum likelihood estimation of abundances. Default: 5000')
    parser.add_option('--junc-alpha', dest='junc_alpha', help='Alpha value for the binomial test used during false positive spliced alignment filtration. Default: 0.001')
    parser.add_option('--small-anchor-fraction', dest='small_anchor_fraction', help='Spliced reads with less than this percent of their length on each side of\
                                                                                      the junction are considered suspicious and are candidates for filtering prior to assembly. Default: 0.09.')
    parser.add_option('--overhang-tolerance', dest='overhang_tolerance', help='The number of bp allowed to enter the intron of a transcript when determining if a \
                                                                                read or another transcript is mappable to/compatible with it. The default is 8 bp based on the default bowtie/TopHat parameters.')
    parser.add_option('--max-bundle-length', dest='max_bundle_length', help='Maximum genomic length of a given bundle" help="Default: 3,500,000bp')
    parser.add_option('--max-bundle-frags', dest='max_bundle_frags', help='Sets the maximum number of fragments a locus may have before being skipped. Skipped loci are listed in skipped.gtf. Default: 1,000,000')
    parser.add_option('--min-intron-length', dest='min_intron_length', help='Minimal allowed intron size. Default: 50')
    parser.add_option('--trim-3-avgcov-thresh', dest='trim_three_avgcov_thresh', help='Minimum average coverage required to attempt 3prime trimming. Default: 10')
    parser.add_option('--trim-3-dropoff-frac', dest='trim_three_dropoff_frac', help='The fraction of average coverage below which to trim the 3prime end of an assembled transcript. Default: 0.1')

    # Bias correction options.
    parser.add_option('-b', dest='do_bias_correction', action="store_true", help='Providing Cufflinks with a multifasta file via this option instructs it to run our new bias detection and correction algorithm which can significantly improve accuracy of transcript abundance estimates.')
    parser.add_option('', '--index', dest='index', help='The path of the reference genome')
    parser.add_option('', '--ref_file', dest='ref_file', help='The reference dataset from the history')

    # Global model (for trackster).
    parser.add_option('', '--global_model', dest='global_model_file', help='Global model used for computing on local data')

    (options, args) = parser.parse_args()

    # If doing bias correction, set/link to sequence file.
    if options.do_bias_correction:
        if options.ref_file:
            # Sequence data from history.
            # Create symbolic link to ref_file so that index will be created in working directory.
            seq_path = "ref.fa"
            os.symlink(options.ref_file, seq_path)
        else:
            if not os.path.exists(options.index):
                stop_err('Reference genome %s not present, request it by reporting this error.' % options.index)
            seq_path = options.index

    # Build command.

    # Base; always use quiet mode to avoid problems with storing log output.
    cmd = "cufflinks -q --no-update-check"

    # Add options.
    if options.max_intron_len:
        cmd += (" -I %i" % int(options.max_intron_len))
    if options.min_isoform_fraction:
        cmd += (" -F %f" % float(options.min_isoform_fraction))
    if options.pre_mrna_fraction:
        cmd += (" -j %f" % float(options.pre_mrna_fraction))
    if options.num_threads:
        cmd += (" -p %i" % int(options.num_threads))
    if options.GTF:
        cmd += (" -G %s" % options.GTF)
    if options.compatible_hits_norm:
        cmd += (" --compatible-hits-norm")
    if options.GTFguide:
        cmd += (" -g %s" % options.GTFguide)
        cmd += (" --3-overhang-tolerance %i" % int(options.three_overhang_tolerance))
        cmd += (" --intron-overhang-tolerance %i" % int(options.intron_overhang_tolerance))
    if options.no_faux_reads:
        cmd += (" --no-faux-reads")
    if options.multi_read_correct:
        cmd += (" -u")

    if options.library_type and options.library_type != 'auto':
        cmd += (" --library-type %s" % options.library_type)
    if options.mask_file:
        cmd += (" --mask-file %s" % options.mask_file)
    if options.inner_mean_dist:
        cmd += (" -m %i" % int(options.inner_mean_dist))
    if options.inner_dist_std_dev:
        cmd += (" -s %i" % int(options.inner_dist_std_dev))
    if options.max_mle_iterations:
        cmd += (" --max-mle-iterations %i" % int(options.max_mle_iterations))
    if options.junc_alpha:
        cmd += (" --junc-alpha %f" % float(options.junc_alpha))
    if options.small_anchor_fraction:
        cmd += (" --small-anchor-fraction %f" % float(options.small_anchor_fraction))
    if options.overhang_tolerance:
        cmd += (" --overhang-tolerance %i" % int(options.overhang_tolerance))
    if options.max_bundle_length:
        cmd += (" --max-bundle-length %i" % int(options.max_bundle_length))
    if options.max_bundle_frags:
        cmd += (" --max-bundle-frags %i" % int(options.max_bundle_frags))
    if options.min_intron_length:
        cmd += (" --min-intron-length %i" % int(options.min_intron_length))
    if options.trim_three_avgcov_thresh:
        cmd += (" --trim-3-avgcov-thresh %i" % int(options.trim_three_avgcov_thresh))
    if options.trim_three_dropoff_frac:
        cmd += (" --trim-3-dropoff-frac %f" % float(options.trim_three_dropoff_frac))

    if options.do_bias_correction:
        cmd += (" -b %s" % seq_path)
    if options.no_effective_length_correction:
        cmd += (" --no-effective-length-correction")
    if options.no_length_correction:
        cmd += (" --no-length-correction")

    # Add input files.
    cmd += " " + options.input

    # Run command and handle output.
    try:
        # Run command.
        with tempfile.NamedTemporaryFile(dir=".") as tmp_stderr:
            returncode = subprocess.call(args=cmd, stderr=tmp_stderr, shell=True)

            # Error checking.
            if returncode != 0:
                # Get stderr, allowing for case where it's very large.
                buffsize = 1048576
                stderr = ''
                with open(tmp_stderr.name) as tmp_stderr2:
                    try:
                        while True:
                            stderr += tmp_stderr2.read(buffsize)
                            if not stderr or len(stderr) % buffsize != 0:
                                break
                    except OverflowError:
                        pass
                raise Exception(stderr)

            # Read standard error to get total map/upper quartile mass.
            total_map_mass = -1
            with open(tmp_stderr.name, 'r') as tmp_stderr2:
                for line in tmp_stderr2:
                    if line.lower().find("map mass") >= 0 or line.lower().find("upper quartile") >= 0:
                        total_map_mass = float(line.split(":")[1].strip())
                        break

        # If there's a global model provided, use model's total map mass
        # to adjust FPKM + confidence intervals.
        if options.global_model_file:
            # Global model is simply total map mass from original run.
            with open(options.global_model_file, 'r') as global_model_file:
                global_model_total_map_mass = float(global_model_file.readline())

            # Ratio of global model's total map mass to original run's map mass is
            # factor used to adjust FPKM.
            fpkm_map_mass_ratio = total_map_mass / global_model_total_map_mass

            # Update FPKM values in transcripts.gtf file.
            with open("transcripts.gtf", 'r') as transcripts_file:
                with tempfile.NamedTemporaryFile(dir=".", delete=False) as new_transcripts_file:
                    for line in transcripts_file:
                        fields = line.split('\t')
                        attrs = parse_gff_attributes(fields[8])
                        attrs["FPKM"] = str(float(attrs["FPKM"]) * fpkm_map_mass_ratio)
                        attrs["conf_lo"] = str(float(attrs["conf_lo"]) * fpkm_map_mass_ratio)
                        attrs["conf_hi"] = str(float(attrs["conf_hi"]) * fpkm_map_mass_ratio)
                        fields[8] = gff_attributes_to_str(attrs, "GTF")
                        new_transcripts_file.write("%s\n" % '\t'.join(fields))
            shutil.move(new_transcripts_file.name, "transcripts.gtf")

        # TODO: update expression files as well.

        # Set outputs. Transcript and gene expression handled by wrapper directives.
        shutil.move("transcripts.gtf", options.assembled_isoforms_output_file)
        if total_map_mass > -1:
            with open("global_model.txt", 'w') as f:
                f.write("%f\n" % total_map_mass)
    except Exception as e:
        stop_err('Error running cufflinks: %s' % e)


if __name__ == "__main__":
    __main__()
