#!/usr/bin/env perl
use warnings;
use strict;

=head1 NAME

windowmasker_to_bed.pl - Convert WindowMasker output to BED format

=head1 SYNOPSIS

  windowmasker_to_bed.pl < ${input.interval} > ${output.bed}

=head1 DESCRIPTION

This scripts converts the interval output from WindowMasker ustat
to bed 4 format.

=head1 VERSION

Last update: 2017-04-28

=cut


#====================
# Libraries
#====================
use Getopt::Long;
use Pod::Usage;
use Carp;

sub main {
  parse_arguments();

  my $chrom;
  my $idx = 0;

  while ( defined(my $line = <>) ) {
    next if ($line =~ /^\s*$/x);

    if ($line =~ /^>(\S+)/x) {
      $chrom = $1;

    } elsif ($line =~ /(\d+)\s+\-\s+(\d+)/x) {
      # interval coordinates are 0-based
      my $start = $1;
      my $end = $2 + 1;
      my $name = sprintf("wm_%d", $idx++);

      print join("\t", $chrom, $start, $end, $name), "\n";

    } else {
      croak("Invalid line: ${line}");
    }
  }

  return;
}

main();



#====================
# Parse command-line arguments
#====================
sub parse_arguments {
  my $help = 0;

  GetOptions('help|?' => \$help) or usage();

  pod2usage({ verbose => 2 }) if ($help);

  return;
}

sub usage {
  my $msg = shift;

  pod2usage({ verbose => 1, message => $msg || "" });
  exit 1;
}
