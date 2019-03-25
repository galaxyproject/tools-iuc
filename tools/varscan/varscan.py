#!/usr/bin/env python3
from __future__ import print_function

import argparse
import io
import os
import subprocess
import sys
import tempfile
import time
from contextlib import ExitStack
from functools import partial
from threading import Thread

import pysam


class VariantCallingError (RuntimeError):
    """Exception class for issues with samtools and varscan subprocesses."""

    def __init__(self, message=None, call='', error=''):
        self.message = message
        self.call = call.strip()
        self.error = error.strip()

    def __str__(self):
        if self.message is None:
            return ''
        if self.error:
            msg_header = '"{0}" failed with:\n{1}\n\n'.format(
                self.call, self.error
            )
        else:
            msg_header = '{0} failed.\n'
            'No further information about this error is available.\n\n'.format(
                self.call
            )
        return msg_header + self.message


class VarScanCaller (object):
    def __init__(self, ref_genome, bam_input_files,
                 max_depth=None,
                 min_mapqual=None, min_basequal=None,
                 threads=1, verbose=False, quiet=True
                 ):
        self.ref_genome = ref_genome
        self.bam_input_files = bam_input_files
        self.max_depth = max_depth
        self.min_mapqual = min_mapqual
        self.min_basequal = min_basequal
        self.threads = threads
        self.verbose = verbose
        self.quiet = quiet

        with pysam.FastaFile(ref_genome) as ref_fa:
            self.ref_contigs = ref_fa.references
            self.ref_lengths = ref_fa.lengths

        self.pileup_engine = ['samtools', 'mpileup']
        self.varcall_engine = ['varscan', 'somatic']
        self.requires_stdout_redirect = False
        self.TemporaryContigVCF = partial(
            tempfile.NamedTemporaryFile,
            mode='wb', suffix='', delete=False, dir=os.getcwd()
        )
        self.tmpfiles = []

    def _get_pysam_pileup_args(self):
        param_dict = {}
        if self.max_depth is not None:
            param_dict['max_depth'] = self.max_depth
        if self.min_mapqual is not None:
            param_dict['min_mapping_quality'] = self.min_mapqual
        if self.min_basequal is not None:
            param_dict['min_base_quality'] = self.min_basequal
        param_dict['compute_baq'] = False
        param_dict['stepper'] = 'samtools'
        return param_dict

    def varcall_parallel(self, normal_purity=None, tumor_purity=None,
                         min_coverage=None,
                         min_var_count=None,
                         min_var_freq=None, min_hom_freq=None,
                         p_value=None, somatic_p_value=None,
                         threads=None, verbose=None, quiet=None
                         ):
        if not threads:
            threads = self.threads
        if verbose is None:
            verbose = self.verbose
        if quiet is None:
            quiet = self.quiet
        # mapping of method parameters to varcall engine command line options
        varcall_engine_option_mapping = [
            ('--normal-purity', normal_purity),
            ('--tumor-purity', tumor_purity),
            ('--min-coverage', min_coverage),
            ('--min-reads2', min_var_count),
            ('--min-var-freq', min_var_freq),
            ('--min-freq-for-hom', min_hom_freq),
            ('--p-value', p_value),
            ('--somatic-p-value', somatic_p_value),
            ('--min-avg-qual', self.min_basequal)
        ]
        varcall_engine_options = []
        for option, value in varcall_engine_option_mapping:
            if value is not None:
                varcall_engine_options += [option, str(value)]
        pileup_engine_options = ['-B']
        if self.max_depth is not None:
            pileup_engine_options += ['-d', str(self.max_depth)]
        if self.min_mapqual is not None:
            pileup_engine_options += ['-q', str(self.min_mapqual)]
        if self.min_basequal is not None:
            pileup_engine_options += ['-Q', str(self.min_basequal)]

        # Create a tuple of calls to samtools mpileup and varscan for
        # each contig. The contig name is stored as the third element of
        # that tuple.
        # The calls are stored in the reverse order of the contig list so
        # that they can be popped off later in the original order
        calls = [(
            self.pileup_engine + pileup_engine_options + [
                '-r', contig + ':',
                '-f', self.ref_genome
            ] + self.bam_input_files,
            self.varcall_engine + [
                '-', '{out}', '--output-vcf', '1', '--mpileup', '1'
            ] + varcall_engine_options,
            contig
        ) for contig in self.ref_contigs[::-1]]

        if verbose:
            print('Starting variant calling ..')

        # launch subprocesses and monitor their status
        subprocesses = []
        error_table = {}
        tmp_io_started = []
        tmp_io_finished = []
        self.tmpfiles = []

        def enqueue_stderr_output(out, stderr_buffer):
            for line in iter(out.readline, b''):
                # Eventually we are going to print the contents of
                # the stderr_buffer to sys.stderr so we can
                # decode things here using its encoding.
                # We do a 'backslashreplace' just to be on the safe side.
                stderr_buffer.write(line.decode(sys.stderr.encoding,
                                                'backslashreplace'))
            out.close()

        try:
            while subprocesses or calls:
                while calls and len(subprocesses) < threads:
                    # There are still calls waiting for execution and we
                    # have unoccupied threads so we launch a new combined
                    # call to samtools mpileup and the variant caller.

                    # pop the call arguments from our call stack
                    call = calls.pop()
                    # get the name of the contig that this call is going
                    # to work on
                    contig = call[2]
                    # Based on the contig name, generate a readable and
                    # file system-compatible prefix and use it to create
                    # a named temporary file, to which the call output
                    # will be redirected.
                    # At the moment we create the output file we add it to
                    # the list of all temporary output files so that we can
                    # remove it eventually during cleanup.
                    call_out = self.TemporaryContigVCF(
                        prefix=''.join(
                            c if c.isalnum() else '_' for c in contig
                        ) + '_',
                    )
                    # maintain a list of variant call outputs
                    # in the order the subprocesses got launched
                    tmp_io_started.append(call_out.name)
                    if self.requires_stdout_redirect:
                        # redirect stdout to the temporary file just created
                        stdout_p2 = call_out
                        stderr_p2 = subprocess.PIPE
                    else:
                        # variant caller wants to write output to file directly
                        stdout_p2 = subprocess.PIPE
                        stderr_p2 = subprocess.STDOUT
                        call[1][call[1].index('{out}')] = call_out.name
                        call_out.close()
                    # for reporting purposes, join the arguments for the
                    # samtools and the variant caller calls into readable
                    # strings
                    c_str = (' '.join(call[0]), ' '.join(call[1]))
                    error_table[c_str] = [io.StringIO(), io.StringIO()]
                    # start the subprocesses
                    p1 = subprocess.Popen(
                        call[0],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE
                    )
                    p2 = subprocess.Popen(
                        call[1],
                        stdin=p1.stdout,
                        stdout=stdout_p2,
                        stderr=stderr_p2
                    )
                    # subprocess bookkeeping
                    subprocesses.append((c_str, p1, p2, call_out, contig))
                    # make sure our newly launched call does not block
                    # because its buffers fill up
                    p1.stdout.close()
                    t1 = Thread(
                        target=enqueue_stderr_output,
                        args=(p1.stderr, error_table[c_str][0])
                    )
                    t2 = Thread(
                        target=enqueue_stderr_output,
                        args=(
                            p2.stderr
                            if self.requires_stdout_redirect else
                            p2.stdout,
                            error_table[c_str][1]
                        )
                    )
                    t1.daemon = t2.daemon = True
                    t1.start()
                    t2.start()

                    if verbose:
                        print(
                            'Calling variants for contig: {0}'
                            .format(call[2])
                        )

                # monitor all running calls to see if any of them are done
                for call, p1, p2, ofo, contig in subprocesses:
                    p1_stat = p1.poll()
                    p2_stat = p2.poll()
                    if p1_stat is not None or p2_stat is not None:
                        # There is an outcome for this process!
                        # Lets see what it is
                        if p1_stat or p2_stat:
                            print()
                            print(
                                error_table[call][0].getvalue(),
                                error_table[call][1].getvalue(),
                                file=sys.stderr
                            )
                            raise VariantCallingError(
                                'Variant Calling for contig {0} failed.'
                                .format(contig),
                                call='{0} | {1}'.format(call[0], call[1])
                            )
                        if p1_stat == 0 and p2_stat is None:
                            # VarScan is not handling the no output from
                            # samtools mpileup situation correctly so maybe
                            # that's the issue here
                            last_words = error_table[call][1].getvalue(
                            ).splitlines()[-4:]
                            if len(last_words) < 4 or any(
                                not msg.startswith('Input stream not ready')
                                for msg in last_words
                            ):
                                break
                                # lets give this process a bit more time
                            # VarScan is waiting for input it will never
                            # get, stop it.
                            p2.terminate()
                            subprocesses.remove((call, p1, p2, ofo, contig))
                            ofo.close()
                            break
                        if p2_stat == 0:
                            # Things went fine.
                            # maintain a list of variant call outputs
                            # that finished successfully (in the order
                            # they finished)
                            tmp_io_finished.append(ofo.name)
                            if verbose:
                                print()
                                print('Contig {0} finished.'.format(contig))
                            if not quiet:
                                print()
                                print(
                                    'stderr output from samtools mpileup/'
                                    'bcftools:'.upper(),
                                    file=sys.stderr
                                )
                                print(
                                    error_table[call][0].getvalue(),
                                    error_table[call][1].getvalue(),
                                    file=sys.stderr
                                )
                            # Discard the collected stderr output from
                            # the call, remove the call from the list of
                            # running calls and close its output file.
                            del error_table[call]
                            subprocesses.remove((call, p1, p2, ofo, contig))
                            # Closing the output file is important or we
                            # may hit the file system limit for open files
                            # if there are lots of contigs.
                            ofo.close()
                            break
                # wait a bit in between monitoring cycles
                time.sleep(2)
        finally:
            for call, p1, p2, ofo, contig in subprocesses:
                # make sure we do not leave running subprocesses behind
                for proc in (p1, p2):
                    try:
                        proc.terminate()
                    except Exception:
                        pass
                # close currently open files
                ofo.close()
            # store the files with finished content in the order that
            # the corresponding jobs were launched
            self.tmpfiles = [f for f in tmp_io_started if f in tmp_io_finished]
            # Make sure remaining buffered stderr output of
            # subprocesses does not get lost.
            # Currently, we don't screen this output for real errors,
            # but simply rewrite everything.
            if not quiet and error_table:
                print()
                print(
                    'stderr output from samtools mpileup/bcftools:'.upper(),
                    file=sys.stderr
                )
                for call, errors in error_table.items():
                    print(' | '.join(call), ':', file=sys.stderr)
                    print('-' * 20, file=sys.stderr)
                    print('samtools mpileup output:', file=sys.stderr)
                    print(errors[0].getvalue(), file=sys.stderr)
                    print('varscan somatic output:', file=sys.stderr)
                    print(errors[1].getvalue(), file=sys.stderr)

    def _add_ref_contigs_to_header(self, header):
        for chrom, length in zip(self.ref_contigs, self.ref_lengths):
            header.add_meta(
                'contig',
                items=[('ID', chrom), ('length', length)]
            )

    def _add_filters_to_header(self, header):
        varscan_fpfilters = {
            'VarCount': 'Fewer than {min_var_count2} variant-supporting reads',
            'VarFreq': 'Variant allele frequency below {min_var_freq2}',
            'VarAvgRL':
                'Average clipped length of variant-supporting reads < '
                '{min_var_len}',
            'VarReadPos': 'Relative average read position < {min_var_readpos}',
            'VarDist3':
                'Average distance to effective 3\' end < {min_var_dist3}',
            'VarMMQS':
                'Average mismatch quality sum for variant reads > '
                '{max_var_mmqs}',
            'VarMapQual':
                'Average mapping quality of variant reads < {min_var_mapqual}',
            'VarBaseQual':
                'Average base quality of variant reads < {min_var_basequal}',
            'Strand':
                'Strand representation of variant reads < {min_strandedness}',
            'RefAvgRL':
                'Average clipped length of ref-supporting reads < '
                '{min_ref_len}',
            'RefReadPos':
                'Relative average read position < {min_ref_readpos}',
            'RefDist3':
                'Average distance to effective 3\' end < {min_ref_dist3}',
            'RefMapQual':
                'Average mapping quality of reference reads < '
                '{min_ref_mapqual}',
            'RefBaseQual':
                'Average base quality of ref-supporting reads < '
                '{min_ref_basequal}',
            'RefMMQS':
                'Average mismatch quality sum for ref-supporting reads > '
                '{max_ref_mmqs}',
            'MMQSdiff':
                'Mismatch quality sum difference (var - ref) > '
                '{max_mmqs_diff}',
            'MinMMQSdiff':
                'Mismatch quality sum difference (var - ref) < '
                '{max_mmqs_diff}',
            'MapQualDiff':
                'Mapping quality difference (ref - var) > {max_mapqual_diff}',
            'MaxBAQdiff':
                'Average base quality difference (ref - var) > '
                '{max_basequal_diff}',
            'ReadLenDiff':
                'Average supporting read length difference (ref - var) > '
                '{max_relative_len_diff}',
        }
        for filter_id, description in varscan_fpfilters.items():
            header.filters.add(filter_id, None, None, description)

    def _add_indel_info_flag_to_header(self, header):
        header.info.add(
            'INDEL', 0, 'Flag', 'Indicates that the variant is an INDEL'
        )

    def _standardize_format_fields(self, header):
        """Add standard FORMAT key declarations to a VCF header."""

        format_field_specs = [
            # pysam will not add quotes around single-word descriptions,
            # which is a violation of the VCF specification.
            # To work around this we are using two words as the
            # GT format description (instead of the commonly used "Genotype").
            # This could be changed once pysam learns to quote correctly.
            ('GT', '1', 'String', 'Genotype code'),
            ('GQ', '1', 'Integer', 'Genotype quality'),
            ('DP', '1', 'Integer', 'Read depth'),
            ('AD', 'R', 'Integer', 'Read depth for each allele'),
            ('ADF', 'R', 'Integer',
             'Read depth for each allele on the forward strand'),
            ('ADR', 'R', 'Integer',
             'Read depth for each allele on the reverse strand')
        ]
        # Existing FORMAT declarations cannot be overwritten.
        # The only viable strategy is to mark them as removed,
        # then merge the header with another one containing the
        # correct declarations. This is what is implemented below.
        temp_header = pysam.VariantHeader()
        for spec in format_field_specs:
            temp_header.formats.add(*spec)
            if spec[0] in header.formats:
                header.formats.remove_header(spec[0])
        header.merge(temp_header)

    def _compile_common_header(self, varcall_template, no_filters=False):
        # read the header produced by VarScan for use as a template
        with pysam.VariantFile(varcall_template, 'r') as original_data:
            varscan_header = original_data.header
        # don't propagate non-standard and redundant FORMAT keys written
        # by VarScan
        varscan_header.formats.remove_header('AD')
        varscan_header.formats.remove_header('FREQ')
        varscan_header.formats.remove_header('RD')
        varscan_header.formats.remove_header('DP4')
        # build a new header containing information not written by VarScan
        common_header = pysam.VariantHeader()
        # copy over sample info from the VarScan template
        for sample in varscan_header.samples:
            common_header.samples.add(sample)
        # add reference information
        common_header.add_meta('reference', value=self.ref_genome)
        # change the source information
        common_header.add_meta('source', value='varscan.py')
        # add contig information
        self._add_ref_contigs_to_header(common_header)
        # declare an INDEL flag for record INFO fields
        self._add_indel_info_flag_to_header(common_header)
        # merge in remaining information from the VarScan template
        common_header.merge(varscan_header)
        # add additional FILTER declarations
        if not no_filters:
            # add filter info
            self._add_filters_to_header(common_header)
        # generate standard FORMAT field declarations
        # including a correct AD declaration to prevent VarScan's
        # non-standard one from getting propagated
        self._standardize_format_fields(common_header)
        return common_header

    def _fix_record_gt_fields(self, record):
        """Migrate non-standard genotype fields to standard ones."""

        # The key issue to consider here is that we need to modify
        # genotype field values on a per-sample basis, but htslib
        # reserves memory for the values of all samples upon the first
        # modification of the field in the variant record.
        # => We need to calculate all values first, then fill each
        # genotype field starting with the sample with the largest value.
        new_gt_fields = {'AD': [], 'ADF': [], 'ADR': []}
        # store the current genotype field contents for each sample
        per_sample_gts = record.samples.values()
        # calculate and store the new genotype field values
        for gt_field in per_sample_gts:
            # generate a standard AD field by combining VarScan's
            # RD and non-standard AD fields
            new_gt_fields['AD'].append((gt_field['RD'], gt_field['AD'][0]))
            # split VarScan's DP4 field into the standard fields
            # ADF and ADR
            new_gt_fields['ADF'].append(
                (int(gt_field['DP4'][0]), int(gt_field['DP4'][2]))
            )
            new_gt_fields['ADR'].append(
                (int(gt_field['DP4'][1]), int(gt_field['DP4'][3]))
            )
        # Modify the record's genotype fields.
        # For each field, modify the sample containing the largest
        # value for the field first.
        # Without this precaution we could trigger:
        # "bcf_update_format: Assertion `!fmt->p_free' failed."
        # in vcf.c of htslib resulting in a core dump.
        for field in sorted(new_gt_fields):
            for new_field, sample_fields in sorted(
                zip(new_gt_fields[field], per_sample_gts),
                key=lambda x: max(x[0]),
                reverse=True
            ):
                sample_fields[field] = new_field
        # remove redundant fields
        # FREQ is easy to calculate from AD
        del record.format['FREQ']
        del record.format['RD']
        del record.format['DP4']

    def get_allele_specific_pileup_column_stats(
        self, allele, pile_column, ref_fetch, use_md=True
    ):
        var_reads_plus = var_reads_minus = 0
        sum_mapping_qualities = 0
        sum_base_qualities = 0
        sum_dist_from_center = 0
        sum_dist_from_3prime = 0
        sum_clipped_length = 0
        sum_unclipped_length = 0
        sum_num_mismatches_as_fraction = 0
        sum_mismatch_qualities = 0

        for p in pile_column.pileups:
            # skip reads that don't support the allele we're looking for
            if p.query_position is None:
                continue
            if p.alignment.query_sequence[p.query_position] != allele:
                continue

            if p.alignment.is_reverse:
                var_reads_minus += 1
            else:
                var_reads_plus += 1
            sum_mapping_qualities += p.alignment.mapping_quality
            sum_base_qualities += p.alignment.query_qualities[p.query_position]
            sum_clipped_length += p.alignment.query_alignment_length
            unclipped_length = p.alignment.query_length
            sum_unclipped_length += unclipped_length
            # The following calculations are all in 1-based coordinates
            # with respect to the physical 5'-end of the read sequence.
            if p.alignment.is_reverse:
                read_base_pos = unclipped_length - p.query_position
                read_first_aln_pos = (
                    unclipped_length - p.alignment.query_alignment_end + 1
                )
                read_last_aln_pos = (
                    unclipped_length - p.alignment.query_alignment_start
                )
            else:
                read_base_pos = p.query_position + 1
                read_first_aln_pos = p.alignment.query_alignment_start + 1
                read_last_aln_pos = p.alignment.query_alignment_end
            # Note: the original bam-readcount algorithm uses the less accurate
            # p.alignment.query_alignment_length / 2 as the read center
            read_center = (read_first_aln_pos + read_last_aln_pos) / 2
            sum_dist_from_center += 1 - abs(
                read_base_pos - read_center
            ) / (read_center - 1)
            # Note: the original bam-readcount algorithm uses a more
            # complicated definition of the 3'-end of a read sequence, which
            # clips off runs of bases with a quality scores of exactly 2.
            sum_dist_from_3prime += (
                read_last_aln_pos - read_base_pos
            ) / (unclipped_length - 1)

            sum_num_mismatches = 0
            if use_md:
                try:
                    # see if the read has an MD tag, in which case pysam can
                    # calculate the reference sequence for us
                    aligned_pairs = p.alignment.get_aligned_pairs(
                        matches_only=True, with_seq=True
                    )
                except ValueError:
                    use_md = False
            if use_md:
                # The ref sequence got calculated from alignment and MD tag.
                for qpos, rpos, ref_base in aligned_pairs:
                    # pysam uses lowercase ref bases to indicate mismatches
                    if (
                        ref_base == 'a' or ref_base == 't' or
                        ref_base == 'g' or ref_base == 'c'
                    ):
                        sum_num_mismatches += 1
                        sum_mismatch_qualities += p.alignment.query_qualities[
                            qpos
                        ]
            else:
                # We need to obtain the aligned portion of the reference
                # sequence.
                aligned_pairs = p.alignment.get_aligned_pairs(
                    matches_only=True
                )
                # note that ref bases can be lowercase
                ref_seq = ref_fetch(
                    aligned_pairs[0][1], aligned_pairs[-1][1] + 1
                ).upper()
                ref_offset = aligned_pairs[0][1]

                for qpos, rpos in p.alignment.get_aligned_pairs(
                    matches_only=True
                ):
                    # see if we have a mismatch to the reference at this
                    # position, but note that
                    # - the query base may be lowercase (SAM/BAM spec)
                    # - an '=', if present in the query seq, indicates a match
                    #   to the reference (SAM/BAM spec)
                    # - there cannot be a mismatch with an N in the reference
                    ref_base = ref_seq[rpos - ref_offset]
                    read_base = p.alignment.query_sequence[qpos].upper()
                    if (
                        read_base != '=' and ref_base != 'N'
                    ) and (
                        ref_base != read_base
                    ):
                        sum_num_mismatches += 1
                        sum_mismatch_qualities += p.alignment.query_qualities[
                            qpos
                        ]
            sum_num_mismatches_as_fraction += (
                sum_num_mismatches / p.alignment.query_alignment_length
            )

        var_reads_total = var_reads_plus + var_reads_minus
        if var_reads_total == 0:
            # No stats without reads!
            return None
        avg_mapping_quality = sum_mapping_qualities / var_reads_total
        avg_basequality = sum_base_qualities / var_reads_total
        avg_pos_as_fraction = sum_dist_from_center / var_reads_total
        avg_num_mismatches_as_fraction = (
            sum_num_mismatches_as_fraction / var_reads_total
        )
        avg_sum_mismatch_qualities = sum_mismatch_qualities / var_reads_total
        avg_clipped_length = sum_clipped_length / var_reads_total
        avg_distance_to_effective_3p_end = (
            sum_dist_from_3prime / var_reads_total
        )

        return (
            avg_mapping_quality,
            avg_basequality,
            var_reads_plus,
            var_reads_minus,
            avg_pos_as_fraction,
            avg_num_mismatches_as_fraction,
            avg_sum_mismatch_qualities,
            avg_clipped_length,
            avg_distance_to_effective_3p_end
        )

    def _postprocess_variant_records(self, invcf,
                                     min_var_count2, min_var_count2_lc,
                                     min_var_freq2, max_somatic_p,
                                     max_somatic_p_depth,
                                     min_ref_readpos, min_var_readpos,
                                     min_ref_dist3, min_var_dist3,
                                     min_ref_len, min_var_len,
                                     max_relative_len_diff,
                                     min_strandedness, min_strand_reads,
                                     min_ref_basequal, min_var_basequal,
                                     max_basequal_diff,
                                     min_ref_mapqual, min_var_mapqual,
                                     max_mapqual_diff,
                                     max_ref_mmqs, max_var_mmqs,
                                     min_mmqs_diff, max_mmqs_diff,
                                     **args):
        # set FILTER field according to Varscan criteria
        # multiple FILTER entries must be separated by semicolons
        # No filters applied should be indicated with MISSING

        # since posterior filters are always applied to just one sample,
        # a better place to store the info is in the FT genotype field:
        # can be PASS, '.' to indicate that filters have not been applied,
        # or a semicolon-separated list of filters that failed
        # unfortunately, gemini does not support this field

        with ExitStack() as io_stack:
            normal_reads, tumor_reads = (
                io_stack.enter_context(
                    pysam.Samfile(fn, 'rb')) for fn in self.bam_input_files
            )
            refseq = io_stack.enter_context(pysam.FastaFile(self.ref_genome))
            pileup_args = self._get_pysam_pileup_args()
            for record in invcf:
                if any(len(allele) > 1 for allele in record.alleles):
                    # skip indel postprocessing for the moment
                    yield record
                    continue
                # get pileup for genomic region affected by this variant
                if record.info['SS'] == '2':
                    # a somatic variant => generate pileup from tumor data
                    pile = tumor_reads.pileup(
                        record.chrom, record.start, record.stop,
                        **pileup_args
                    )
                    sample_of_interest = 'TUMOR'
                elif record.info['SS'] in ['1', '3']:
                    # a germline or LOH variant => pileup from normal data
                    pile = normal_reads.pileup(
                        record.chrom, record.start, record.stop,
                        **pileup_args
                    )
                    sample_of_interest = 'NORMAL'
                else:
                    # TO DO: figure out if there is anything interesting to do
                    # for SS status codes 0 (reference) and 5 (unknown)
                    yield record
                    continue

                # apply false-positive filtering a la varscan fpfilter
                # find the variant site in the pileup columns
                for pile_column in pile:
                    if pile_column.reference_pos == record.start:
                        break
                # extract required information
                # overall read depth at the site
                read_depth = pile_column.get_num_aligned()
                assert read_depth > 0
                # no multiallelic sites in varscan
                assert len(record.alleles) == 2
                if record.samples[sample_of_interest]['RD'] > 0:
                    ref_stats, alt_stats = [
                        self.get_allele_specific_pileup_column_stats(
                            allele,
                            pile_column,
                            partial(
                                pysam.FastaFile.fetch, refseq, record.chrom
                            )
                        )
                        for allele in record.alleles
                    ]
                else:
                    ref_stats = None
                    alt_stats = self.get_allele_specific_pileup_column_stats(
                        record.alleles[1],
                        pile_column,
                        partial(pysam.FastaFile.fetch, refseq, record.chrom)
                    )
                ref_count = 0
                if ref_stats:
                    ref_count = ref_stats[2] + ref_stats[3]
                    if ref_stats[1] < min_ref_basequal:
                        record.filter.add('RefBaseQual')
                    if ref_count >= 2:
                        if ref_stats[0] < min_ref_mapqual:
                            record.filter.add('RefMapQual')
                        if ref_stats[4] < min_ref_readpos:
                            record.filter.add('RefReadPos')
                        # ref_stats[5] (avg_num_mismatches_as_fraction
                        # is not a filter criterion in VarScan fpfilter
                        if ref_stats[6] > max_ref_mmqs:
                            record.filter.add('RefMMQS')
                        if ref_stats[7] < min_ref_len:
                            # VarScan fpfilter does not apply this filter
                            # for indels, but there is no reason
                            # not to do it.
                            record.filter.add('RefAvgRL')
                        if ref_stats[8] < min_ref_dist3:
                            record.filter.add('RefDist3')
                if alt_stats:
                    alt_count = alt_stats[2] + alt_stats[3]
                    if (
                        alt_count < min_var_count2_lc
                    ) or (
                        read_depth >= max_somatic_p_depth and
                        alt_count < min_var_count2
                    ):
                        record.filter.add('VarCount')
                    if alt_count / read_depth < min_var_freq2:
                        record.filter.add('VarFreq')
                    if alt_stats[1] < min_var_basequal:
                        record.filter.add('VarBaseQual')
                    if alt_count > min_strand_reads:
                        if (
                            alt_stats[2] / alt_count < min_strandedness
                        ) or (
                            alt_stats[3] / alt_count < min_strandedness
                        ):
                            record.filter.add('Strand')
                    if alt_stats[2] + alt_stats[3] >= 2:
                        if alt_stats[0] < min_var_mapqual:
                            record.filter.add('VarMapQual')
                        if alt_stats[4] < min_var_readpos:
                            record.filter.add('VarReadPos')
                        # alt_stats[5] (avg_num_mismatches_as_fraction
                        # is not a filter criterion in VarScan fpfilter
                        if alt_stats[6] > max_var_mmqs:
                            record.filter.add('VarMMQS')
                        if alt_stats[7] < min_var_len:
                            # VarScan fpfilter does not apply this filter
                            # for indels, but there is no reason
                            # not to do it.
                            record.filter.add('VarAvgRL')
                        if alt_stats[8] < min_var_dist3:
                            record.filter.add('VarDist3')
                    if ref_count >= 2 and alt_count >= 2:
                        if (ref_stats[0] - alt_stats[0]) > max_mapqual_diff:
                            record.filter.add('MapQualDiff')
                        if (ref_stats[1] - alt_stats[1]) > max_basequal_diff:
                            record.filter.add('MaxBAQdiff')
                        mmqs_diff = alt_stats[6] - ref_stats[6]
                        if mmqs_diff < min_mmqs_diff:
                            record.filter.add('MinMMQSdiff')
                        if mmqs_diff > max_mmqs_diff:
                            record.filter.add('MMQSdiff')
                        if (
                            1 - alt_stats[7] / ref_stats[7]
                        ) > max_relative_len_diff:
                            record.filter.add('ReadLenDiff')
                else:
                    # No variant-supporting reads for this record!
                    # This can happen in rare cases because of
                    # samtools mpileup issues, but indicates a
                    # rather unreliable variant call.
                    record.filter.add('VarCount')
                    record.filter.add('VarFreq')
                yield record

    def _indel_flagged_records(self, vcf):
        for record in vcf:
            record.info['INDEL'] = True
            yield record

    def _merge_generator(self, vcf1, vcf2):
        try:
            record1 = next(vcf1)
        except StopIteration:
            for record2 in vcf2:
                yield record2
            return
        try:
            record2 = next(vcf2)
        except StopIteration:
            yield record1
            for record1 in vcf1:
                yield record1
            return
        while True:
            if (record1.start, record1.stop) < (record2.start, record2.stop):
                yield record1
                try:
                    record1 = next(vcf1)
                except StopIteration:
                    yield record2
                    for record2 in vcf2:
                        yield record2
                    return
            else:
                yield record2
                try:
                    record2 = next(vcf2)
                except StopIteration:
                    yield record1
                    for record1 in vcf1:
                        yield record1
                    return

    def merge_and_postprocess(self, snps_out, indels_out=None,
                              no_filters=False, **filter_args):
        temporary_data = self.tmpfiles
        self.tmpfiles = []
        temporary_snp_files = [f + '.snp.vcf' for f in temporary_data]
        temporary_indel_files = [f + '.indel.vcf' for f in temporary_data]

        for f in temporary_data:
            try:
                os.remove(f)
            except Exception:
                pass

        def noop_gen(data, **kwargs):
            for d in data:
                yield d

        if no_filters:
            apply_filters = noop_gen
        else:
            apply_filters = self._postprocess_variant_records

        output_header = self._compile_common_header(
            temporary_snp_files[0],
            no_filters
        )
        if indels_out is None:
            with open(snps_out, 'w') as o:
                o.write(str(output_header).format(**filter_args))
                for snp_f, indel_f in zip(
                    temporary_snp_files, temporary_indel_files
                ):
                    with pysam.VariantFile(snp_f, 'r') as snp_invcf:
                        # fix the input header on the fly
                        # to avoid Warnings from htslib about missing contig
                        # info
                        self._add_ref_contigs_to_header(snp_invcf.header)
                        self._add_filters_to_header(snp_invcf.header)
                        self._add_indel_info_flag_to_header(snp_invcf.header)
                        self._standardize_format_fields(snp_invcf.header)
                        with pysam.VariantFile(indel_f, 'r') as indel_invcf:
                            # fix the input header on the fly
                            # to avoid Warnings from htslib about missing
                            # contig info
                            self._add_ref_contigs_to_header(indel_invcf.header)
                            self._add_filters_to_header(indel_invcf.header)
                            self._add_indel_info_flag_to_header(
                                indel_invcf.header
                            )
                            self._standardize_format_fields(indel_invcf.header)
                            for record in apply_filters(
                                self._merge_generator(
                                    snp_invcf,
                                    self._indel_flagged_records(indel_invcf)
                                ),
                                **filter_args
                            ):
                                self._fix_record_gt_fields(record)
                                o.write(str(record))
                    try:
                        os.remove(snp_f)
                    except Exception:
                        pass
                    try:
                        os.remove(indel_f)
                    except Exception:
                        pass

        else:
            with open(snps_out, 'w') as o:
                o.write(str(output_header).format(**filter_args))
                for f in temporary_snp_files:
                    with pysam.VariantFile(f, 'r') as invcf:
                        # fix the input header on the fly
                        # to avoid Warnings from htslib about missing
                        # contig info and errors because of undeclared
                        # filters
                        self._add_ref_contigs_to_header(invcf.header)
                        self._add_filters_to_header(invcf.header)
                        self._standardize_format_fields(invcf.header)
                        for record in apply_filters(
                            invcf, **filter_args
                        ):
                            self._fix_record_gt_fields(record)
                            o.write(str(record))
                    try:
                        os.remove(f)
                    except Exception:
                        pass
            with open(indels_out, 'w') as o:
                o.write(str(output_header))
                for f in temporary_indel_files:
                    with pysam.VariantFile(f, 'r') as invcf:
                        # fix the input header on the fly
                        # to avoid Warnings from htslib about missing
                        # contig info and errors because of undeclared
                        # filters
                        self._add_ref_contigs_to_header(invcf.header)
                        self._add_filters_to_header(invcf.header)
                        self._add_indel_info_flag_to_header(invcf.header)
                        self._standardize_format_fields(invcf.header)
                        for record in apply_filters(
                            self._indel_flagged_records(invcf), **filter_args
                        ):
                            self._fix_record_gt_fields(record)
                            o.write(str(record))
                    try:
                        os.remove(f)
                    except Exception:
                        pass


def varscan_call(ref_genome, normal, tumor, output_path, **args):
    """Preparse arguments and orchestrate calling and postprocessing."""

    if args.pop('split_output'):
        if '%T' in output_path:
            out = (
                output_path.replace('%T', 'snp'),
                output_path.replace('%T', 'indel')
            )
        else:
            out = (
                output_path + '.snp',
                output_path + '.indel'
            )
    else:
        out = (output_path, None)

    instance_args = {
        k: args.pop(k) for k in [
            'max_depth',
            'min_mapqual',
            'min_basequal',
            'threads',
            'verbose',
            'quiet'
        ]
    }
    varscan_somatic_args = {
        k: args.pop(k) for k in [
            'normal_purity',
            'tumor_purity',
            'min_coverage',
            'min_var_count',
            'min_var_freq',
            'min_hom_freq',
            'somatic_p_value',
            'p_value'
        ]
    }

    v = VarScanCaller(ref_genome, [normal, tumor], **instance_args)
    v.varcall_parallel(**varscan_somatic_args)
    v.merge_and_postprocess(*out, **args)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument(
        'ref_genome',
        metavar='reference_genome',
        help='the reference genome (in fasta format)'
    )
    p.add_argument(
        '--normal',
        metavar='BAM_file', required=True,
        help='the BAM input file of aligned reads from the normal sample'
    )
    p.add_argument(
        '--tumor',
        metavar='BAM_file', required=True,
        help='the BAM input file of aligned reads from the tumor sample'
    )
    p.add_argument(
        '-o', '--ofile', required=True,
        metavar='OFILE', dest='output_path',
        help='Name of the variant output file. With --split-output, the name '
             'may use the %%T replacement token or will be used as the '
             'basename for the two output files to be generated (see '
             '-s|--split-output below).'
    )
    p.add_argument(
        '-s', '--split-output',
        dest='split_output', action='store_true', default=False,
        help='indicate that separate output files for SNPs and indels '
             'should be generated (original VarScan behavior). If specified, '
             '%%T in the --ofile file name will be replaced with "snp" and '
             '"indel" to generate the names of the SNP and indel output '
             'files, respectively. If %%T is not found in the file name, it '
             'will get interpreted as a basename to which ".snp"/".indel" '
             'will be appended.'
    )
    p.add_argument(
        '-t', '--threads',
        type=int, default=1,
        help='level of parallelism'
    )
    p.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='be verbose about progress'
    )
    p.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='suppress output from wrapped tools'
    )
    call_group = p.add_argument_group('Variant calling parameters')
    call_group.add_argument(
        '--normal-purity',
        dest='normal_purity', type=float,
        default=1.0,
        help='Estimated purity of the normal sample (default: 1.0)'
    )
    call_group.add_argument(
        '--tumor-purity',
        dest='tumor_purity', type=float,
        default=1.0,
        help='Estimated purity of the tumor sample (default: 1.0)'
    )
    call_group.add_argument(
        '--max-pileup-depth',
        dest='max_depth', type=int, default=8000,
        help='Maximum depth of generated pileups (samtools mpileup -d option; '
             'default: 8000)'
    )
    call_group.add_argument(
        '--min-basequal',
        dest='min_basequal', type=int,
        default=13,
        help='Minimum base quality at the variant position to use a read '
             '(default: 13)'
    )
    call_group.add_argument(
        '--min-mapqual',
        dest='min_mapqual', type=int,
        default=0,
        help='Minimum mapping quality required to use a read '
             '(default: 0)'
    )
    call_group.add_argument(
        '--min-coverage',
        dest='min_coverage', type=int,
        default=8,
        help='Minimum site coverage required in the normal and in the tumor '
             'sample to call a variant (default: 8)'
    )
    call_group.add_argument(
        '--min-var-count',
        dest='min_var_count', type=int,
        default=2,
        help='Minimum number of variant-supporting reads required to call a '
             'variant (default: 2)'
    )
    call_group.add_argument(
        '--min-var-freq',
        dest='min_var_freq', type=float,
        default=0.1,
        help='Minimum variant allele frequency for calling (default: 0.1)'
    )
    call_group.add_argument(
        '--min-hom-freq',
        dest='min_hom_freq', type=float,
        default=0.75,
        help='Minimum variant allele frequency for homozygous call '
             '(default: 0.75)'
    )
    call_group.add_argument(
        '--p-value',
        dest='p_value', type=float,
        default=0.99,
        help='P-value threshold for heterozygous call (default: 0.99)'
    )
    call_group.add_argument(
        '--somatic-p-value',
        dest='somatic_p_value', type=float,
        default=0.05,
        help='P-value threshold for somatic call (default: 0.05)'
    )
    filter_group = p.add_argument_group('Posterior variant filter parameters')
    filter_group.add_argument(
        '--no-filters',
        dest='no_filters', action='store_true',
        help='Disable all posterior variant filters. '
             'If specified, all following options will be ignored'
    )
    filter_group.add_argument(
        '--min-var-count2',
        dest='min_var_count2', type=int,
        default=4,
        help='Minimum number of variant-supporting reads (default: 4)'
    )
    filter_group.add_argument(
        '--min-var-count2-lc',
        dest='min_var_count2_lc', type=int,
        default=2,
        help='Minimum number of variant-supporting reads when depth below '
             '--somatic-p-depth (default: 2)'
    )
    filter_group.add_argument(
        '--min-var-freq2',
        dest='min_var_freq2', type=float,
        default=0.05,
        help='Minimum variant allele frequency (default: 0.05)'
    )
    filter_group.add_argument(
        '--max-somatic-p',
        dest='max_somatic_p', type=float,
        default=0.05,
        help='Maximum somatic p-value (default: 0.05)'
    )
    filter_group.add_argument(
        '--max-somatic-p-depth',
        dest='max_somatic_p_depth', type=int,
        default=10,
        help='Depth required to run --max-somatic-p filter (default: 10)'
    )
    filter_group.add_argument(
        '--min-ref-readpos',
        dest='min_ref_readpos', type=float,
        default=0.1,
        help='Minimum average relative distance of site from the ends of '
             'ref-supporting reads (default: 0.1)'
    )
    filter_group.add_argument(
        '--min-var-readpos',
        dest='min_var_readpos', type=float,
        default=0.1,
        help='Minimum average relative distance of site from the ends of '
             'variant-supporting reads (default: 0.1)'
    )
    filter_group.add_argument(
        '--min-ref-dist3',
        dest='min_ref_dist3', type=float,
        default=0.1,
        help='Minimum average relative distance of site from the effective '
             '3\'end of ref-supporting reads (default: 0.1)'
    )
    filter_group.add_argument(
        '--min-var-dist3',
        dest='min_var_dist3', type=float,
        default=0.1,
        help='Minimum average relative distance of site from the effective '
             '3\'end of variant-supporting reads (default: 0.1)'
    )
    filter_group.add_argument(
        '--min-ref-len',
        dest='min_ref_len', type=int,
        default=90,
        help='Minimum average trimmed length of reads supporting the ref '
             'allele (default: 90)'
    )
    filter_group.add_argument(
        '--min-var-len',
        dest='min_var_len', type=int,
        default=90,
        help='Minimum average trimmed length of reads supporting the variant '
             'allele (default: 90)'
    )
    filter_group.add_argument(
        '--max-len-diff',
        dest='max_relative_len_diff', type=float,
        default=0.25,
        help='Maximum average relative read length difference (ref - var; '
             'default: 0.25)'
    )
    filter_group.add_argument(
        '--min-strandedness',
        dest='min_strandedness', type=float,
        default=0.01,
        help='Minimum fraction of variant reads from each strand '
             '(default: 0.01)'
    )
    filter_group.add_argument(
        '--min-strand-reads',
        dest='min_strand_reads', type=int,
        default=5,
        help='Minimum allele depth required to run --min-strandedness filter '
             '(default: 5)'
    )
    filter_group.add_argument(
        '--min-ref-basequal',
        dest='min_ref_basequal', type=int,
        default=15,
        help='Minimum average base quality for the ref allele (default: 15)'
    )
    filter_group.add_argument(
        '--min-var-basequal',
        dest='min_var_basequal', type=int,
        default=15,
        help='Minimum average base quality for the variant allele '
             '(default: 15)'
    )
    filter_group.add_argument(
        '--max-basequal-diff',
        dest='max_basequal_diff', type=int,
        default=50,
        help='Maximum average base quality diff (ref - var; default: 50)'
    )
    filter_group.add_argument(
        '--min-ref-mapqual',
        dest='min_ref_mapqual', type=int,
        default=15,
        help='Minimum average mapping quality of reads supporting the ref '
             'allele (default: 15)'
    )
    filter_group.add_argument(
        '--min-var-mapqual',
        dest='min_var_mapqual', type=int,
        default=15,
        help='Minimum average mapping quality of reads supporting the variant '
             'allele (default: 15)'
    )
    filter_group.add_argument(
        '--max-mapqual-diff',
        dest='max_mapqual_diff', type=int,
        default=50,
        help='Maximum average mapping quality difference (ref - var; '
             'default: 50)'
    )
    filter_group.add_argument(
        '--max-ref-mmqs',
        dest='max_ref_mmqs', type=int,
        default=100,
        help='Maximum mismatch quality sum of reads supporting the ref '
             'allele (default: 100)'
    )
    filter_group.add_argument(
        '--max-var-mmqs',
        dest='max_var_mmqs', type=int,
        default=100,
        help='Maximum mismatch quality sum of reads supporting the variant '
             'allele (default: 100)'
    )
    filter_group.add_argument(
        '--min-mmqs-diff',
        dest='min_mmqs_diff', type=int,
        default=0,
        help='Minimum mismatch quality sum difference (var - ref; default: 0)'
    )
    filter_group.add_argument(
        '--max-mmqs-diff',
        dest='max_mmqs_diff', type=int,
        default=50,
        help='Maximum mismatch quality sum difference (var - ref; default: 50)'
    )
    args = vars(p.parse_args())
    varscan_call(**args)
