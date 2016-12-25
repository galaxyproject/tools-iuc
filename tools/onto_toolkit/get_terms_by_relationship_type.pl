#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
   if 0; # not running under some shell

use Carp;
use strict;
use warnings;
use Encode;

use OBO::Parser::OBOParser;

my $my_parser     = OBO::Parser::OBOParser->new();
my $ontology      = $my_parser->work(shift(@ARGV));
my $rel_type_name = shift(@ARGV);
my @relationships = @{$ontology->get_relationships()};

foreach my $r (@relationships) {
	if($r->type() eq $rel_type_name){
		print $r->tail()->id(), "\t";
		print $r->type(),"\t";
		print $r->head()->id(), "\n";
	}
}

exit 0;

__END__
