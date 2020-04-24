#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use lib dirname(__FILE__);
use multialign qw(windows);

##G_windows.pl
#sam   $sam
#out   $dir
#win   $win_num
#help  $help

my  ($proc, $uni, $sam, $ma,$win_num, $dir_uni, $dir_non, $dir, $zip, $help, %sens, %antisens);
$win_num = 1000;

GetOptions (
"sam:s"  =>  \$sam,
"dir=s"  =>  \$dir,
"ma=i"  =>  \$ma,
"uni" => \$uni,
"zip" => \$zip,
"p:i" => \$proc,
"win:i" => \$win_num,
"help"   =>  \$help
);

die "--sam STRING --dir STRING  --win INT [1000] (O for auto)  --uni --ma --zip \n" if $help;

die "No input file specified\n" unless defined($sam);
die "No output specified\n" unless defined($dir);
print STDERR "-----------------------------\n";
print STDERR "creation of Gviz:\nsam_file: $sam\nout_dir: $dir\nwindows number: $win_num\n";

mkdir $dir;
$dir = $dir.'/' unless $dir =~ /(.*)\/$/;

$win_num = 'auto' if $win_num == 0;

if ($uni)
{
    $dir_uni = $dir.'unique/';
    #print STDERR ("DEBUG_G_windows_dir_uni:".$dir_uni."\n");
    mkdir $dir_uni;
    $dir_non = $dir.'non_unique/';
    #print STDERR ("DEBUG_G_windows_dir_non:".$dir_non."\n");
    mkdir $dir_non;  
    windows($sam, $win_num, $dir_non , 1, $dir_uni,$proc, $ma);
}
else { windows($sam, $win_num, $dir , 0, "NA", $proc, $ma);}

if ($zip)
{
    my $zipfile = $1 if $dir =~ /(.*)\/$/;
    #print STDERR ("bash -c 'zip -r $zipfile $dir'"."\n");
    system('zip -r $zipfile $dir');
    if ($uni)
    {
        unlink glob "$dir_uni*";
        unlink glob "$dir_non*";
        rmdir $dir_uni;
        rmdir $dir_non;
    }
    unlink glob "$dir*"; 
    system('rm -r $dir');
}
