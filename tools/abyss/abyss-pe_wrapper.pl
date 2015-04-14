#!/usr/bin/env perl

use strict;
use warnings;
use File::Copy;

die("ERROR: Expected at least 8 arguments; got: @ARGV\n") unless @ARGV >= 8;

# GET ARGS
my $kmer_size=shift @ARGV;
my $min_num_pairs=shift @ARGV;
my $outdir=shift @ARGV;
my $outfile=shift @ARGV;
my $contigs_outfile=shift @ARGV;
my $sam_outfile=shift @ARGV;
my $hist_outfile=shift @ARGV;

# ALL FILES GO IN THIS extra_files_path
unless (-d $outdir) {
    mkdir $outdir or die("Unable to make dir, $outdir\n");
}
chdir $outdir;

# RUN COMMAND
`abyss-pe k=$kmer_size n=$min_num_pairs in='@ARGV' name=abyss 2> $outdir/abyss.stderr > $outdir/abyss.stdout`;
if ($? != 0) {
    unless ( -s "$outdir/abyss-3.hist") { print STDERR "NO CONTIGS WERE PRODUCED!\n" }
    open(IN, "<$outdir/abyss.stdout") or die($!);
    while (<IN>) { print STDERR $_ }
    close IN;
    die("ABORTING\n");
}

# FILTER HISTOGRAM
open (IN, "<$outdir/abyss.stderr") or die($!);
open (OUT, ">$outfile") or die($!);
while (my $line=<IN>) {
    my @chars=split(//, $line);
    my $filter=0;
    foreach my $char (@chars) {
        if (ord($char) >= 129) {
            $filter=1;
            last;
        }
    }
    print OUT $line unless $filter;
}
close IN;
close OUT;
unlink("$outdir/abyss.stderr");

# OUTPUT INFO LINE TEXT
open(IN, "<$outdir/abyss.stdout") or die($!);
while (<IN>) {
    if (/^Assembled \d+ k\-mer in \d+ contigs/) {
        print;
        last;
    }
}
close IN;

# GALAXY DOESN'T WANT GZIPPED DATAFILES
run("gunzip $outdir/abyss-3.sam.gz");

# MOVE OUTFILES
mv_outfile("$outdir/abyss-contigs.fa", $contigs_outfile);
mv_outfile("$outdir/abyss-3.sam", $sam_outfile);
mv_outfile("$outdir/coverage.hist", $hist_outfile);
exit;

sub run {
    my $cmd=shift;
    my $output=`$cmd 2>&1`;
    if ($? != 0) {
        print STDERR "ERROR RUNNING COMMAND: $cmd\n";
        die($output);
    }
    return $output;
}


sub mv_outfile {
    my ($src,$dest)=@_;
    # if dest defined and src exist, then move outfiles to galaxy-specified location
    if ( $dest ne 'None' and -f $src ) {
        unlink($dest);
        move($src,$dest);
    }
}
__END__
