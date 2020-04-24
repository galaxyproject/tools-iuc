#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use multialign qw(align histogram rpms_rpkm);

my (@files, $fni, $fastq, $dir, $raw, $rawU, $zip, $fin,$frej, $ref,$proc, $fout, $funi, $sam, $dup, $mis, $help, %genome_hits, %duplicates, $repartition, $png_rep, $rep_size,$aligner, $noD);
$mis = 3;
my ($ma, $pi, $mi, $bo, $min, $max) =(0,0,0,0,0,0);

GetOptions (
"fin=s"  =>  \$fin ,
"funi" => \$funi,
"min=i" =>  \$min ,
"max=i" =>  \$max ,
"ref=s"  =>  \$ref ,
"dir:s" => \$dir,
"zip" => \$zip,
"dup" => \$dup,
"ma=i" => \$ma,
"pi=i" => \$pi,
"mi=i" => \$mi,
"bo=i" => \$bo,
"rep" => \$repartition,
"noD" => \$noD,
"mismatches=i" => \$mis,
"help"   =>  \$help,
"aligner:s" => \$aligner,
"p=i" => \$proc,
);

die "--fin STRING --ref STRING --dir STRING  --aligner Bowtie2|BWA --mis INT [3] --ma INT --pi INT[0] --p --mi INT[0] --bo INT[0] --funi --dup --zip\n" if $help;
die "No input file specified\n" unless defined($fin);
die "No reference specified\n" unless defined($ref);
die "No aligner specified\n" unless defined($aligner);
chomp($aligner,$fin,$ref);
if (defined ($dir))
{
    my @suffix = ('.fastq','.ref', '.fq', '.dat', '.fa','.fas','.fasta', '.txt');
    my ($name,$path,$suffix) = fileparse($fin,@suffix);
    my ($name_r,$path_r,$suffix_r) = fileparse($ref,@suffix);
    $dir = $dir.'/' unless $dir =~ /(.*)\/$/;
    mkdir $dir;
    $fout = $dir.$name.'_'.$name_r.'.fastq';  push(@files, $fout);
    if ($funi) {$fni = $dir.$name.'_unique_'.$name_r.'.fastq';   push(@files, $fni);}
    else {$fni = '/dev/null';}
    $frej = $dir.$name.'_rejected_'.$name_r.'.fastq';   push(@files, $frej);
    $sam = $dir.$name.'_'.$name_r.'_sorted.sam';   push(@files, $sam);
    $raw = $dir.$name.'_'.$name_r.'_reads_counts.txt';   push(@files, $raw);
    $rawU = $dir.$name.'_'.$name_r.'_unique_reads_counts.txt';   push(@files, $rawU);
    
    if ($repartition)
    {
        $repartition = $dir.$name.'_'.$name_r.'_distribution.txt';  push(@files, $repartition);
        $png_rep = $dir.$name.'_'.$name_r.'_distribution.png';
        $rep_size = $dir.$name.'_'.$name_r.'_size_distribution/' unless $noD;
    }
}
else{  die "No output directory specified\n";}

print STDERR "-----------------------------\n";
print STDERR "Align and count:\nfastq_in: $fin\nreference: $ref\nmismatches: $mis\n";

#alignement to the first reference
my @return = align($ref, $mis, $fin, $fout, $fni, $frej, $sam, \%duplicates, \%genome_hits, $aligner, $proc);
my $uniqueNum  = $return[4];
my $num = $return[3];
my $size_mappedHashR = $return[-2];
my $size_mapped_spe_HashR = $return[-1];
print $num."\t".$uniqueNum."\n";
#print number of duplicates and hits number

if ($repartition)
{
    my ($rpm, $percentage, $total) =(0,0,0);
    $total += $_ foreach values %{$size_mappedHashR};
    open (my $rep, '>'.$repartition) or die "Cannot open: $! \n";
    print $rep "size\tnumber\tpercentage";
    print $rep "\trpm" if $ma != 0;
    print $rep "\n";
    if ($min != 0 && $max!= 0)
    {
    for my $i ($min..$max) { ${$size_mappedHashR}{$i} = 0 unless exists(${$size_mappedHashR}{$i});}
    }
    if ($min != 0 && $max!= 0)
    {
    for my $i ($min..$max) { ${$size_mappedHashR}{$i} = 0 unless exists(${$size_mappedHashR}{$i});}
    }
    foreach my $k (sort{$a <=> $b} keys (%{$size_mappedHashR}))
    {
        $percentage = $size_mappedHashR->{$k} / $total * 100;
        $rpm = $size_mappedHashR->{$k} / $ma * 1000000 if $ma != 0;
        print $rep "$k\t$size_mappedHashR->{$k}\t";
        printf $rep "%.2f",$percentage;
        print $rep "\t$rpm";
        print $rep "\n";
    }

    eval 
    {
        histogram($size_mappedHashR, $png_rep, $total); 
        print STDERR ("INFO eval histo size:".$size_mappedHashR."\n");
        print STDERR ("INFO eval histo size:".$png_rep."\n");
        print STDERR ("INFO eval histo size:".$total."\n");
    };

    unless ($noD)
    {
        mkdir $rep_size;
        foreach my $k (keys %{$size_mapped_spe_HashR})
        {
            $k =~ tr/:/-/;
            $repartition = $rep_size.$k.'.txt';
            $png_rep = $rep_size.$k.'.png';
            my $hashRef = ${$size_mapped_spe_HashR}{$k};
            open ($rep, '>'.$repartition) or die "Cannot open $repartition: $!\n";
            $total = 0;
            $total += $_ foreach values %{$hashRef};
            print $rep "size\tnumber\tpercentage";
            print $rep "\trpm" if $ma != 0;
            print $rep "\n";
            if ($min != 0 && $max!= 0)
            {
                for my $i ($min..$max) { ${$hashRef}{$i} = 0 unless exists(${$hashRef}{$i});}
            }

            foreach my $l (sort{$a <=> $b} keys (%{$hashRef}))
            {
                $percentage = 0;
                $percentage =  $hashRef->{$l} / $total * 100 unless $total == 0;
                $rpm = $hashRef->{$l} / $ma * 1000000 if $ma != 0;
                print $rep "$l\t$hashRef->{$l}\t";
                printf $rep "%.2f",$percentage;
                print $rep "\t$rpm" if $ma != 0;
                print $rep "\n";
            }
            eval { histogram($size_mapped_spe_HashR->{$k}, $png_rep, $total); };
        }
    }
}

$ma = $return[3] if $ma == 0;
rpms_rpkm($return[0], $return[2], $ma, $raw, $pi, $mi, $bo);
rpms_rpkm($return[1], $return[2] ,$ma,  $rawU, $pi, $mi, $bo) if ($funi);

#print number of duplicates and hits number
if (defined $dup)
{
    $dup = $dir.'dup_mapnum.txt';  push(@files, $dup);
    open(my $tab,">".$dup) or die "Cannot open output txt file: $!\n";
    print $tab "sequence\tcount\tmapnum\n";
    foreach my $k (sort {$duplicates{$b} <=> $duplicates{$a}}keys %duplicates)
    {
        $duplicates{$k} = 0 unless exists($duplicates{$k});
        $genome_hits{$k} = 0 unless exists($genome_hits{$k});
        print $tab $k."\t".$duplicates{$k}."\t".$genome_hits{$k}."\n";
    }
    close $tab;
}

if ($zip)
{
    my $zipfile = $1 if $dir =~ /(.*)\/$/;
    system('zip -r $zipfile $dir');
    unlink @files; rmdir $dir;
    system('rm -r $rep_size unless $noD');
}
