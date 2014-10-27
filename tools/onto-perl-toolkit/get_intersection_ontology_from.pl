#!/usr/bin/perl -w

use OBO::Parser::OBOParser;
use OBO::Util::Ontolome;


use Carp;
use strict;
use warnings;

my @input_files = @ARGV;

my $parse_ont = OBO::Parser::OBOParser->new();
my $ome = OBO::Util::Ontolome->new();

my @input_ontos;
foreach(@input_files)
{
	my $onto = $parse_ont->work($_) or warn "ontology not defined, check input files\n";
	push @input_ontos, $onto;
}

my $intersection_ontology = $ome->intersection(@input_ontos);
$intersection_ontology->export(\*STDOUT, "obo");

exit 0; 

