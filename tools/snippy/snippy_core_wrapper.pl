#!/usr/bin/env perl

#--------------------------------------
#
#        snippy_core_wrapper.pl
#
# This is an intermediary script between snippy-core.xml and snippy-core
# It:
#   - Copys the supplied zipped snippy output files to the working dir
#   - Untars them to their datafile name
#   - Builds the snippy-core command
#   - Runs the snippy-core command
#
#--------------------------------------

use warnings;
use strict;
use File::Copy;
use File::Basename;

my(@Options, $indirs, $noref);
setOptions();

my @list_of_dirs = split /\s+/, $indirs;

#The list of final directories to be passed to snippy-core will be stored here.
my @snippy_outs;

foreach my $d (@list_of_dirs){
  print STDERR "$d\n";
  my $bn = basename($d);
  my ($name, $dir, $ext) = fileparse($d, '\..*');
  copy($d, $bn);
  print STDERR "$d, $bn, $name, $dir, $ext\n";
  `tar -xf $bn`;
}

my $test_list = `ls -d */`;
$test_list =~ s/\///g;
print STDERR "$test_list\n";

@snippy_outs = split /\s+/, $test_list;


my $commandline = "snippy-core ";

$commandline .= "--noref " if $noref;
$commandline .= join(" ", @snippy_outs);
print STDERR "snippy-core commandline: $commandline\n";

my $ok = system($commandline);

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"noref!", VAR=>\$noref, DEFAULT=>0, DESC=>"Don't include the reference in the alignment."},
    {OPT=>"indirs=s", VAR=>\$indirs, DEFAULT=>"", DESC=>"A whitespace delimited list of the snippy output zipped dirs."},
  );

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] -i inputfile > ... \n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
