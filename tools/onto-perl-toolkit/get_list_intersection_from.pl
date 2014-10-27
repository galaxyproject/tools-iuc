#!/usr/bin/perl -w

use Carp;
use strict;
use warnings;

# Takes in first argument
my ($input_file1, $input_file2 ) = @_;

# Opens the input files
	open my $fh, "<", $ARGV[0] or die "could not open $ARGV[0]: $!"; 
	my @list1 = <$fh>; 


	open my $fh1, "<", $ARGV[1] or die "could not open $ARGV[1]: $!";
	my @list2 = <$fh1>;

# generates an intersection of terms from the first and second arguments
my @intersection_terms;
foreach (@list1)
{
	my $itemlist1 = $_;
	
	foreach (@list2)
	{
		my $itemlist2 = $_;
		
		if ($itemlist1 eq $itemlist2)
		{
			push @intersection_terms, $itemlist2;
		}
	}	
}
print @intersection_terms, "\n";


exit 0;
