#!/usr/local/bin/env perl
use strict;
use warnings;

my $cmd_string = join (" ",@ARGV);
my $results = `$cmd_string`;
my @files = split("\n",$results);
foreach my $thisLine (@files)
{
	if ($thisLine =~ /Created /)
	{
		$thisLine =~ /[\w|\.]+$/;
		$thisLine =$&;
		print "outfile: $thisLine\n";
	}
	else
	{
		print $thisLine,"\n";
	}
}
