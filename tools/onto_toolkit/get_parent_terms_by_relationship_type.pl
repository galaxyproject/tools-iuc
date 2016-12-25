#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
   if 0; # not running under some shell


use Carp;
use strict;
use warnings;

use OBO::Parser::OBOParser;

my $my_parser = OBO::Parser::OBOParser->new();
my $ontology  = $my_parser->work(shift(@ARGV));
my $term_id   = shift(@ARGV);
my $rel_id    = shift(@ARGV);

my $target_term = $ontology->get_term_by_id($term_id);
my @rels        = @{$ontology->get_relationships_by_source_term($target_term, $rel_id)};

foreach my $r (@rels) {
      # print "rel: ", $r->id(), "\n";
       print $r->head()->id();
      # print "tail: ", $r->tail()->id(), "\n\n";
}

exit 0;
__END__






