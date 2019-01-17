#!/usr/bin/env perl

# fasta-stats
# written by torsten.seemann@monash.edu 
# oct 2012

use strict;
use warnings;
use List::Util qw(sum min max);

# stat storage

my $n=0;
my $seq = '';
my %stat;
my @len;

# MAIN LOOP collecting sequences

while (my $line = <ARGV>) {
  chomp $line;
  if ($line =~ m/^\s*>/) {
    process($seq) if $n;
    $n++;
    $seq='';
  }
  else {
    $seq .= $line;
  }  
}

process($seq) if $n;

# sort length array 
# (should use hash here for efficiency with huge no of short reads?)

@len = sort { $a <=> $b } @len;

# compute more stats

$stat{'num_seq'} = scalar(@len);

if (@len) {
  $stat{'num_bp'} = sum(@len);
  $stat{'len_min'} = $len[0];
  $stat{'len_max'} = $len[-1];
  $stat{'len_median'} = $len[int(@len/2)];
  $stat{'len_mean'} = int( $stat{'num_bp'} / $stat{'num_seq'} ); 
  # calculate n50

  $stat{'len_N50'} = 0;
  my $cum=0;
  my $thresh = int 0.5 * $stat{'num_bp'};
  for my $i (0 .. $#len) {
    $cum += $len[$i];
    if ($cum >= $thresh) {
      $stat{'len_N50'} = $len[$i];
      last;
    }
  }
}

#calculate GC content
$stat{'num_bp_not_N'} = $stat{'num_G'} + $stat{'num_C'} + $stat{'num_A'} + $stat{'num_T'};
$stat{'GC_content'} = ($stat{'num_G'} + $stat{'num_C'}) / $stat{'num_bp_not_N'}*100;

# print stats as .tsv

for my $name (sort keys %stat) {
    if ($name =~ m/GC_content/){
        printf "%s\t%0.1f\n", $name, $stat{$name};
    } else {
        printf "%s\t%s\n", $name, $stat{$name};
    }
}

# run for each sequence

sub process {
  my($s) = @_;
  # base composition
  for my $x (qw(A G T C N)) {
    my $count = $s =~ s/$x/$x/gi;
    $stat{"num_$x"} += $count;
  }
  # keep list of all lengths encountered
  push @len, length($s);    
}

