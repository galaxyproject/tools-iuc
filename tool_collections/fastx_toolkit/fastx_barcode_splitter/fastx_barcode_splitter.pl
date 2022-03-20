#!/usr/bin/env perl

#    FASTX-toolkit - FASTA/FASTQ preprocessing tools.
#    Copyright (C) 2009-2013  A. Gordon (assafgordon@gmail.com)
#
#   Lance Parsons (lparsons@princeton.edu)
#   3/21/2011 - Modified to accept separate index file for barcodes
#   4/6/2011 - Modified to cleanup bad barcode identifiers (esp. useful for Galaxy)
#   4/28/2016 - Modified summary output to remove file paths and add comment
#               character '#'

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Affero General Public License as
#   published by the Free Software Foundation, either version 3 of the
#   License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use IO::Handle;
use Data::Dumper;
use Getopt::Long;
use Carp;

##
## This program splits a FASTQ/FASTA file into several smaller files,
## Based on barcode matching.
##
## run with "--help" for usage information
##
## Assaf Gordon <assafgordon@gmail.com> , 11sep2008

# Forward declarations
sub load_barcode_file ($);
sub parse_command_line ;
sub match_sequences ;
sub mismatch_count($$) ;
sub print_results;
sub open_and_detect_input_format;
sub open_index_and_detect_input_format($);
sub read_index_record;
sub read_record;
sub write_record($);
sub usage();

# Global flags and arguments,
# Set by command line argumens
my $barcode_file ;
my $barcodes_at_eol = 0 ;
my $barcodes_at_bol = 0 ;
my $index_read_file ;
my $exact_match = 0 ;
my $allow_partial_overlap = 0;
my $allowed_mismatches = 1;
my $newfile_suffix = '';
my $newfile_prefix  ;
my $quiet = 0 ;
my $debug = 0 ;
my $fastq_format = 1;
my $index_fastq_format = 1;
my $read_id_check_strip_characters = 1;

# Global variables
# Populated by 'create_output_files'
my %filenames;
my %files;
my %counts = ( 'unmatched' => 0 );
my $barcodes_length;
my @barcodes;
my $input_file_io;


# The Four lines per record in FASTQ format.
# (when using FASTA format, only the first two are used)
my $seq_name;
my $seq_bases;
my $seq_name2;
my $seq_qualities;

# Values used for index read file
my $index_seq_name;
my $index_seq_bases;
my $index_seq_name2;
my $index_seq_qualities;

#
# Start of Program
#
parse_command_line ;

load_barcode_file ( $barcode_file ) ;

open_and_detect_input_format;

if (defined $index_read_file) {open_index_and_detect_input_format ( $index_read_file );}

match_sequences ;

print_results unless $quiet;

#
# End of program
#

sub parse_command_line {
    my $help;

    usage() if (scalar @ARGV==0);

    my $result = GetOptions ( "bcfile=s" => \$barcode_file,
                  "eol"  => \$barcodes_at_eol,
                  "bol"  => \$barcodes_at_bol,
                  "idxfile=s"  => \$index_read_file,
                  "idxidstrip=i" => \$read_id_check_strip_characters,
                  "exact" => \$exact_match,
                  "prefix=s" => \$newfile_prefix,
                  "suffix=s" => \$newfile_suffix,
                  "quiet" => \$quiet,
                  "partial=i" => \$allow_partial_overlap,
                  "debug" => \$debug,
                  "mismatches=i" => \$allowed_mismatches,
                  "help" => \$help
    ) ;

    usage() if ($help);

    die "Error: barcode file not specified (use '--bcfile [FILENAME]')\n" unless defined $barcode_file;
    die "Error: prefix path/filename not specified (use '--prefix [PATH]')\n" unless defined $newfile_prefix;

    if (! defined $index_read_file) {
        if ($barcodes_at_bol == $barcodes_at_eol) {
            die "Error: can't specify both --eol & --bol\n" if $barcodes_at_eol;
            die "Error: must specify either --eol or --bol or --idxfile\n" ;
        }
    }
    elsif ($barcodes_at_bol || $barcodes_at_eol) {
        die "Error: Must specify only one of --idxfile, --eol, or --bol";
    }

    die "Error: invalid for value partial matches (valid values are 0 or greater)\n" if $allow_partial_overlap<0;

    $allowed_mismatches = 0 if $exact_match;

    die "Error: invalid value for mismatches (valid values are 0 or more)\n" if ($allowed_mismatches<0);

    die "Error: partial overlap value ($allow_partial_overlap) bigger than " .
        "max. allowed mismatches ($allowed_mismatches)\n" if ($allow_partial_overlap > $allowed_mismatches);


    exit unless $result;
}



#
# Read the barcode file
#
sub load_barcode_file ($) {
    my $filename = shift or croak "Missing barcode file name";

    open BCFILE,"<$filename" or die "Error: failed to open barcode file ($filename)\n";
    while (<BCFILE>) {
        next if m/^#/;
        chomp;
        my ($ident, $barcode) = split('\t') ;

        $barcode = uc($barcode);

        # Sanity checks on the barcodes
        die "Error: bad data at barcode file ($filename) line $.\n" unless defined $barcode;
        die "Error: bad barcode value ($barcode) at barcode file ($filename) line $.\n"
            unless $barcode =~ m/^[AGCT]+$/;

        # Cleanup Identifiers (only allow alphanumeric, replace others with dash '-')
        $ident =~ s/[^A-Za-z0-9]/-/g;
        die "Error: bad identifier value ($ident) at barcode file ($filename) line $. (must be alphanumeric)\n"
            unless $ident =~ m/^\w+$/;

        die "Error: badcode($ident, $barcode) is shorter or equal to maximum number of " .
            "mismatches ($allowed_mismatches). This makes no sense. Specify fewer mismatches.\n"
                if length($barcode)<=$allowed_mismatches;

        $barcodes_length = length($barcode) unless defined $barcodes_length;
        die "Error: found barcodes in different lengths. this feature is not supported yet.\n"
            unless $barcodes_length == length($barcode);

        push @barcodes, [$ident, $barcode];

        if ($allow_partial_overlap>0) {
            foreach my $i (1 .. $allow_partial_overlap) {
                substr $barcode, ($barcodes_at_bol)?0:-1, 1, '';
                push @barcodes, [$ident, $barcode];
            }
        }
    }
    close BCFILE;

    if ($debug) {
        print STDERR "barcode\tsequence\n";
        foreach my $barcoderef (@barcodes) {
            my ($ident, $seq) = @{$barcoderef};
            print STDERR $ident,"\t", $seq ,"\n";
        }
    }
}

# Create one output file for each barcode.
# (Also create a file for the dummy 'unmatched' barcode)
sub create_output_files {
    my %barcodes = map { $_->[0] => 1 } @barcodes; #generate a uniq list of barcode identifiers;
    $barcodes{'unmatched'} = 1 ;

    foreach my $ident (keys %barcodes) {
        my $new_filename = $newfile_prefix . $ident . $newfile_suffix;
        $filenames{$ident} = $new_filename;
        open my $file, ">$new_filename" or die "Error: failed to create output file ($new_filename)\n";
        $files{$ident} = $file ;
    }
}

sub match_sequences {

    my %barcodes = map { $_->[0] => 1 } @barcodes; #generate a uniq list of barcode identifiers;
    $barcodes{'unmatched'} = 1 ;

    #reset counters
    foreach my $ident ( keys %barcodes ) {
        $counts{$ident} = 0;
    }

    create_output_files;

    # Read file FASTQ file
    # split accotding to barcodes
    while ( read_record ) {
        chomp $seq_name;
        chomp $seq_bases;
        if (defined $index_read_file) {
            read_index_record() or die "Error: Unable to read index sequence for sequence name ($seq_name), check to make sure the file lengths match.\n";
            chomp $index_seq_name;
            chomp $index_seq_bases;

            # Assert that the read ids match
            my $seq_name_match = &strip_read_id($seq_name);
            my $index_seq_name_match = &strip_read_id($index_seq_name);
            if ($seq_name_match ne $index_seq_name_match) {
                die "Error: Index sequence name ($index_seq_name) does not match sequence name ($seq_name)\n";
            }

        }

        print STDERR "sequence $seq_bases: \n" if $debug;

        my $best_barcode_mismatches_count = $barcodes_length;
        my $best_barcode_ident = undef;

        #Try all barcodes, find the one with the lowest mismatch count
        foreach my $barcoderef (@barcodes) {
            my ($ident, $barcode) = @{$barcoderef};

            # Get DNA fragment (in the length of the barcodes)
            # The barcode will be tested only against this fragment
            # (no point in testing the barcode against the whole sequence)
            my $sequence_fragment;
            if ($barcodes_at_bol) {
                $sequence_fragment = substr $seq_bases, 0, $barcodes_length;
            } elsif ($barcodes_at_eol) {
                $sequence_fragment = substr $seq_bases, - $barcodes_length;
            } else {
                $sequence_fragment = substr $index_seq_bases, 0, $barcodes_length;
            }

            my $mm = mismatch_count($sequence_fragment, $barcode) ;

            # if this is a partial match, add the non-overlap as a mismatch
            # (partial barcodes are shorter than the length of the original barcodes)
            $mm += ($barcodes_length - length($barcode));

            if ( $mm < $best_barcode_mismatches_count ) {
                $best_barcode_mismatches_count = $mm ;
                $best_barcode_ident = $ident ;
            }
        }

        $best_barcode_ident = 'unmatched'
            if ( (!defined $best_barcode_ident) || $best_barcode_mismatches_count>$allowed_mismatches) ;

        print STDERR "sequence $seq_bases matched barcode: $best_barcode_ident\n" if $debug;

        $counts{$best_barcode_ident}++;

        #get the file associated with the matched barcode.
        #(note: there's also a file associated with 'unmatched' barcode)
        my $file = $files{$best_barcode_ident};

        write_record($file);
    }
}

# Strip end of readids when matching to avoid mismatch between read 1, 2, 3, etc.
sub strip_read_id {
    my $read_id = shift;
    my $stripped_read_id = $read_id;
    if ($read_id_check_strip_characters) {
        if ($read_id =~ /@([^:]+):([0-9]+):([^:]+):([0-9]+):([0-9]+):([0-9]+):([0-9]+) ([0-9]+):([YN]):([0-9]+):([ACGT]+){0,1}/) { # CASAVA 1.8+
            my @parts = split(/ /,$read_id);
            $stripped_read_id = $parts[0];
        } else { # CASAVA 1.7 and earlier
            $stripped_read_id = substr($read_id, 0, length($read_id)-$read_id_check_strip_characters);
        }
    }
    return $stripped_read_id;
}

#Quickly calculate hamming distance between two strings
#
#NOTE: Strings must be same length.
#      returns number of different characters.
#see  http://www.perlmonks.org/?node_id=500235
sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }

sub print_results
{
    print "# Barcode\tCount\n";
    my $total = 0 ;
    foreach my $ident (sort keys %counts) {
        print $ident, "\t", $counts{$ident},"\n";
        $total += $counts{$ident};
    }
    print "total\t",$total,"\n";
}

sub read_record
{
    $seq_name = $input_file_io->getline();

    return undef unless defined $seq_name; # End of file?

    $seq_bases = $input_file_io->getline();
    die "Error: bad input file, expecting line with sequences\n" unless defined $seq_bases;

    # If using FASTQ format, read two more lines
    if ($fastq_format) {
        $seq_name2 = $input_file_io->getline();
        die "Error: bad input file, expecting line with sequence name2\n" unless defined $seq_name2;

        $seq_qualities = $input_file_io->getline();
        die "Error: bad input file, expecting line with quality scores\n" unless defined $seq_qualities;
    }
    return 1;
}

sub write_record($)
{
    my $file = shift;

    croak "Bad file handle" unless defined $file;

    print $file $seq_name,"\n";
    print $file $seq_bases,"\n";

    #if using FASTQ format, write two more lines
    if ($fastq_format) {
        print $file $seq_name2;
        print $file $seq_qualities;
    }
}

sub open_and_detect_input_format
{
    $input_file_io = new IO::Handle;
    die "Failed to open STDIN " unless $input_file_io->fdopen(fileno(STDIN),"r");

    # Get the first characeter, and push it back
    my $first_char = $input_file_io->getc();
    $input_file_io->ungetc(ord $first_char);

    if ($first_char eq '>') {
        # FASTA format
        $fastq_format = 0 ;
        print STDERR "Detected FASTA format\n" if $debug;
    } elsif ($first_char eq '@') {
        # FASTQ format
        $fastq_format = 1;
        print STDERR "Detected FASTQ format\n" if $debug;
    } else {
        die "Error: unknown file format. First character = '$first_char' (expecting > or \@)\n";
    }
}

sub open_index_and_detect_input_format($) {
    my $filename = shift or croak "Missing index read file name";

    open IDXFILE,"<$filename" or die "Error: failed to open index read file ($filename)\n";

    # Get the first line, and reset file pointer
    my $first_line = <IDXFILE>;
    my $first_char = substr($first_line, 0, 1);
    seek(IDXFILE, 0, 0);

    if ($first_char eq '>') {
        # FASTA format
        $index_fastq_format = 0 ;
        print STDERR "Detected FASTA format for index file\n" if $debug;
    } elsif ($first_char eq '@') {
        # FASTQ format
        $index_fastq_format = 1;
        print STDERR "Detected FASTQ format for index file\n" if $debug;
    } else {
        die "Error: unknown index file format. First character = '$first_char' (expecting > or \@)\n";
    }
}

sub read_index_record
{
    $index_seq_name = <IDXFILE>;

    return undef unless defined $index_seq_name; # End of file?

    $index_seq_bases = <IDXFILE>;
    die "Error: bad input file, expecting line with sequences\n" unless defined $index_seq_bases;

    # If using FASTQ format, read two more lines
    if ($index_fastq_format) {
        $index_seq_name2    = <IDXFILE>;
        die "Error: bad input file, expecting line with sequence name2\n" unless defined $index_seq_name2;

        $index_seq_qualities = <IDXFILE>;
        die "Error: bad input file, expecting line with quality scores\n" unless defined $index_seq_qualities;
    }
    return 1;
}

sub usage()
{
    print<<EOF;
Barcode Splitter, by Assaf Gordon (gordon\@cshl.edu), 11sep2008

This program reads FASTA/FASTQ file and splits it into several smaller files,
Based on barcode matching.
FASTA/FASTQ data is read from STDIN (format is auto-detected.)
Output files will be writen to disk.
Summary will be printed to STDOUT.

usage: $0 --bcfile FILE --prefix PREFIX [--suffix SUFFIX] [--bol|--eol|--idxfile]
         [--mismatches N] [--exact] [--partial N] [--idxidstrip N]
         [--help] [--quiet] [--debug]

Arguments:

--bcfile FILE  - Barcodes file name. (see explanation below.)
--prefix PREFIX  - File prefix. will be added to the output files. Can be used
      to specify output directories.
--suffix SUFFIX  - File suffix (optional). Can be used to specify file
      extensions.
--bol    - Try to match barcodes at the BEGINNING of sequences.
      (What biologists would call the 5' end, and programmers
      would call index 0.)
--eol    - Try to match barcodes at the END of sequences.
      (What biologists would call the 3' end, and programmers
      would call the end of the string.)
--idxfile FILE  - Read barcodes from separate index file (fasta or fastq)
      NOTE: one of --bol, --eol, --idxfile must be specified,
           but not more than one.
--idxidstrip N  - When using index file, strip this number of characters
      from the end of the sequence id before matching.
      Automatically detects CASAVA 1.8 format and strips at a
      space in the id, use 0 to disable this.
      (Default is 1).
--mismatches N  - Max. number of mismatches allowed. default is 1.
--exact    - Same as '--mismatches 0'. If both --exact and --mismatches
      are specified, '--exact' takes precedence.
--partial N  - Allow partial overlap of barcodes. (see explanation below.)
      (Default is not partial matching)
--quiet    - Don't print counts and summary at the end of the run.
      (Default is to print.)
--debug    - Print lots of useless debug information to STDERR.
--help    - This helpful help screen.

Example (Assuming 's_2_100.txt' is a FASTQ file, 'mybarcodes.txt' is
the barcodes file):

  \$ cat s_2_100.txt | $0 --bcfile mybarcodes.txt --bol --mismatches 2 \\
  --prefix /tmp/bla_ --suffix ".txt"

Barcode file format
-------------------
Barcode files are simple text files. Each line should contain an identifier
(descriptive name for the barcode), and the barcode itself (A/C/G/T),
separated by a TAB character. Example:

    #This line is a comment (starts with a 'number' sign)
    BC1 GATCT
    BC2 ATCGT
    BC3 GTGAT
    BC4 TGTCT

For each barcode, a new FASTQ file will be created (with the barcode's
identifier as part of the file name). Sequences matching the barcode
will be stored in the appropriate file.

Running the above example (assuming "mybarcodes.txt" contains the above
barcodes), will create the following files:
  /tmp/bla_BC1.txt
  /tmp/bla_BC2.txt
  /tmp/bla_BC3.txt
  /tmp/bla_BC4.txt
  /tmp/bla_unmatched.txt
The 'unmatched' file will contain all sequences that didn't match any barcode.

Barcode matching
----------------

** Without partial matching:

Count mismatches between the FASTA/Q sequences and the barcodes.
The barcode which matched with the lowest mismatches count (providing the
count is small or equal to '--mismatches N') 'gets' the sequences.

Example (using the above barcodes):
Input Sequence:
GATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG

Matching with '--bol --mismatches 1':
GATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG
GATCT (1 mismatch, BC1)
ATCGT (4 mismatches, BC2)
GTGAT (3 mismatches, BC3)
TGTCT (3 mismatches, BC4)

This sequence will be classified as 'BC1' (it has the lowest mismatch count).
If '--exact' or '--mismatches 0' were specified, this sequence would be
classified as 'unmatched' (because, although BC1 had the lowest mismatch count,
it is above the maximum allowed mismatches).

Matching with '--eol' (end of line) does the same, but from the other side
of the sequence.

** With partial matching (very similar to indels):

Same as above, with the following addition: barcodes are also checked for
partial overlap (number of allowed non-overlapping bases is '--partial N').

Example:
Input sequence is ATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG
(Same as above, but note the missing 'G' at the beginning.)

Matching (without partial overlapping) against BC1 yields 4 mismatches:
ATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG
GATCT (4 mismatches)

Partial overlapping would also try the following match:
-ATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG
GATCT (1 mismatch)

Note: scoring counts a missing base as a mismatch, so the final
mismatch count is 2 (1 'real' mismatch, 1 'missing base' mismatch).
If running with '--mismatches 2' (meaning allowing upto 2 mismatches) - this
seqeunce will be classified as BC1.

EOF

exit 1;
}
