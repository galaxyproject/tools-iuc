#!/usr/bin/env perl
#

use strict;
use warnings;
use Getopt::Long;

my($otuFile, $mFile, $outputFile, $notaxidFile, $label, $krona, $help, @fds, @new_tax);
my $taxid = 0;

GetOptions( "otuTable=s"  => \$otuFile,
            "query=s"     => \$mFile,
            "outfile=s"   => \$outputFile,
            "krona=s"     => \$krona,
            "label=s"     => \$label,
            "h|help"      => \$help,
            "taxid!"      => \$taxid,
            "notaxidfile:s" => \$notaxidFile,    ) or die "Unknown option.\n";

if($help){
    help();
}

#No help has been requested.  Sanity checks....
if(!$label){
    $label = "Unspecified";
}


#Check that the OTU table is defined and present.
if(!defined($otuFile)){
    die "The OTU table is not defined\n";
}elsif(!-e $otuFile){
    die "The OUT table, $otuFile, does not exist\n";
}

#Check that the query file is defined and present.
if(!defined($mFile)){
    die "The mapseq output file (query) file is not defined\n";
}elsif(!-e $mFile){
    die "The mapseq file, $mFile, does not exist\n";
}

#Check that the output file is defined.
if(!defined($outputFile)){
    die "The output file file is not defined\n";
}

#This gives us an idea about what the mapseq file format.

# mapseq v1.0 (Nov 13 2016)
#querydbhitbitscoreidentitymatchesmismatchesgapsquery_startquery_enddbhit_startdbhit_endSILVA
#contig--1213616/899-1048AY664005.1.12231410.98666668148110150170319D_0__Bacteria;D_1__Cyanobacteria;D_2__Cyanobacteria;D_3__SubsectionI;D_4__FamilyI;D_5__Synechococcus;D_6__uncultured Synechococcus sp.
#contig--1126892/723-858EU805193.1.12931340.99264705135100136821957D_0__Bacteria
#contig--243158/4-340JQ611080.1.9103250.99107143332113370335D_0__Archaea;D_1__Euryarchaeota;D_2__Thermoplasmata;D_3__Thermoplasmatales;D_4__Marine Group II;D_5__uncultured archaeon;D_6__uncultured archaeon
#contig--1205795/1003-7EU802400.1.14959810.99395779876049970993D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Oceanospirillales;D_4__SAR86 clade
#contig--1213616/623-456EU802790.1.12661410.94674557160631168146314D_0__Bacteria;D_1__Cyanobacteria;D_2__Cyanobacteria;D_3__SubsectionI;D_4__FamilyI
#contig--2489731/672-501LURT01000153.4866.63211410.9506173115471016212951456D_0__Archaea;D_1__Euryarchaeota;D_2__Thermoplasmata;D_3__Thermoplasmatales;D_4__Marine Group II
#contig--151250/871-2AACY020187844.576.19958430.9850402585613018700869D_0__Archaea;D_1__Euryarchaeota;D_2__Thermoplasmata;D_3__Thermoplasmatales;D_4__Marine Group II
#contig--1206060/3-639EF574940.1.14506080.998360636091006108401450D_0__Bacteria;D_1__Cyanobacteria;D_2__Cyanobacteria;D_3__SubsectionI;D_4__FamilyI

#Need to count up when we see a taxonomy string multiple times.
open(M, "<", $mFile) or die "Could not open mapseq results file $mFile:[$!]\n";
my $taxCount;
my $tax="";
while(<M>){
    if(/^#/){
    #Skip counts
	next;
    }else{
    #Pull out the fields that we need
	chomp;
	my @line = split(/\t/);
	if(!$line[13]){
	    $tax = "Unclassified";
	} else {
	    $tax= $line[13];
	    until ($tax !~/\_\_$/) {
		@fds=split (/\;/, $tax);
		splice @fds, -1;
		$tax= join(";",@fds);
	    }
	}
	$taxCount->{$tax}->{count}++;
    }
}
close(M) or die "Could not close filehandle on mapseq file\n.";



#Now we do not need to read the whole OTU table, just those that we have found
#in mapseq results.
#my $otuFile = "consensus_taxonomy_7_levels.otu";
open(O, "<", $otuFile) or die "Could not open OTU file $otuFile:[$!]\n";
while(<O>){
    chomp;
  #Simple two column table, OTU code and taxonomy string.
    my ($otu, $tax, $taxid) = split(/\t/, $_, 3);
  #print "$otu | $tax";
  #Have we seen this tax string? Store the OTU
    if (defined($taxCount->{$tax})){
	$taxCount->{$tax}->{otu} = $otu;
	if ($taxid){
	    $taxCount->{$tax}->{taxid} = $taxid;
	}
    }
}
close(O) or die "Failed to close filehande on OTU table.\n";




#now check all taxonomy strings have an otu. If this fails,
#something is very wrong.
foreach my $tax (keys %{$taxCount}){
    if(!defined($taxCount->{$tax}->{otu})){
	die "Fatal, |$tax| has not got an OTU code assigned\n";
    }
}

if ($taxid){
    if ($label =~ /^UNITE/) {
	open(T, ">", $outputFile) or die "Could not open 'taxid_file':[$!]\n";
	print T "# Constructed from biom file\n# OTU ID\t".$label."\ttaxonomy\n";
    } else{
	open(T, ">", $outputFile) or die "Could not open 'taxid_file':[$!]\n";
	print T "# Constructed from biom file\n# OTU ID\t".$label."\ttaxonomy\ttaxid\n";
    }
    open(R, ">", $notaxidFile ) or die "Could not open $outputFile:[$!]\n";
    print R "# Constructed from biom file\n# OTU ID\t".$label."\ttaxonomy\n";
}else{
    open(R, ">", $outputFile ) or die "Could not open $outputFile:[$!]\n";
    print R "# Constructed from biom file\n# OTU ID\t".$label."\ttaxonomy\n";
}

#Print the file out
foreach my $tax (sort{$a cmp $b} keys %{$taxCount}){
    if ($taxid){
	print T $taxCount->{$tax}->{otu}."\t".sprintf("%.1f", $taxCount->{$tax}->{count})."\t".$tax."\t".$taxCount->{$tax}->{taxid}."\n";
	print R $taxCount->{$tax}->{otu}."\t".sprintf("%.1f", $taxCount->{$tax}->{count})."\t".$tax."\n";
    }else{
	print R $taxCount->{$tax}->{otu}."\t".sprintf("%.1f", $taxCount->{$tax}->{count})."\t".$tax."\n";
    }
}

close(R) or die "Failed to close open filehandle on outputfile\n";
close(T) or die "Failed to close open filehandle on notaxidfile\n";

if($krona){
    open(K, ">", $krona ) or die "Could not open $krona:[$!]\n";
    foreach my $tax (sort{$a cmp $b} keys %{$taxCount}){
	my $taxMod = $tax;
	$taxMod =~ s/\D\_\d{1}\_\_//g;
	my @tax = split(/\;/, $taxMod);
	$taxMod = join("\t", @tax);
	print K $taxCount->{$tax}->{count}."\t".$taxMod."\n";
    }
    close(K) or die "Failed to close open filehandle on outputfile\n";
}



sub help {

    print<<EOF;

$0 --otuTable onsensus_taxonomy_7_levels.otu --query ERR1234567.outfile --outfile ERR1234567.tsv

Options
-------
  otuTable : <string>, the OTU table produced for the taxonomies found in the reference databases that was used with MAPseq.
  query    : <string>, the output from the MAPseq that assigns a taxonomy to a sequence.
  outfile  : <string>, the file storing the tsv file.
  krona    : <string>, output file name for the Krona text. (Optional).
  label    : <string>, lable to add to the top of the outfile OTU table.
  h|help   : prints this message.

EOF

exit 1;
}
