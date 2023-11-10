#!/usr/bin/perl
use strict;
use warnings;

# Written by Gerrit Timmerhaus (gerrit.timmerhaus@biologie.uni-freiburg.de).
# Changes included by Kristian Ullrich, Per Wilhelmsson and Romy Petroll.

# Script to extract all detected domains out of a hmmsearch results file and classify the families of all used proteins based on these domains.
# The classification depends on a table which contains all known classification rules for the protein families of interest and on specific coverage values defined for every domain.
# The script provides three outputs, namely output.1, output.2 and output.3. The output files are tables in ";"-delimited format.
# The structure of output.1 is: "sequence ID ; TAP family ; number of classifications ; domains". 
# Output.3 shares in principle the same structure as output.1, except that subfamilies are considered. ("sequence ID ; TAP family ; Subfamily ; number of classifications ; domains")
# The superior TAP family is specified first, followed by the subfamily. If a TAP family has no subfamily, the TAP family is specified first and then a "-". 
# The structure of output.2 is: "TAP family";"number of detected proteins".
# More than one entry for a protein is possible because the classification rules may allow more than one classification.
#
# The script must be startet with the arguments <hmmsearch output file> <classification rules> <output classifications file> <output family statistics file> <output subfamily classifications file> <"filter" if desired>

if (!@ARGV or ($ARGV [0] eq "-h") or ($ARGV [0] eq "-help")) {
	print "Usage: extract.and.classify.pl <hmmsearch output file> <classification rules> <output classifications file> <output family statistics file> <output subfamily classifications file> <\"filter\" (if desired)>\n\n";
	exit;
}

# hmmsearch_output: domtblout file
my $hmmsearch_output = $ARGV [0];
# decision_table: rules file
my $decision_table = $ARGV [1];
# family_classifications: output.1
my $family_classifications = $ARGV [2]; 
# family_statistics: output.2
my $family_statistics = $ARGV [3];
# subfamily_classifications: output.3
my $subfamily_classifications = $ARGV [4];
# domspec_cuts: coverage values file
my $domspec_cuts = $ARGV [5];
# gene_model_filte: filter for ARATH and ORYSA
my $gene_model_filter = $ARGV [6];

if ($family_statistics eq "") {
	print "Usage: extract.and.classify.pl <hmmsearch output file> <classification rules> <output classifications file> <output family statistics file> <output subfamil classifications file> <\"filter\" (if desired)>\n\n";
	exit;
}

if ($gene_model_filter and $gene_model_filter eq "filter") {
	print "\nGene model filter is activated. It only works for TAIR (Arabidopsis) and TIGR (Rice) proteins up to now\n";
}

# Array where the $hmmsearch-output/domtblout file will be stored
my @output = ();
# Array with domain-specific coverage values
my @cuts = ();
# Array with rules
my @dec_table = ();
# Counter for the number of detected domains in the hmmsearch output file
my $entry_counter = 0;
# Containes the actual result for a query sequence
my $akt_entry = "";
# Used to define query entry to ignore similar domains
my $whole_entry = ""; 
# Includes the final entries after ignoring similar domains
my @results_of_extraction = ();
# Used to define query entry to ignore similar domains
my $extracted_domain = "";
# Used to define query entry to ignore similar domains
my $present = "";
# Used to define query entry to ignore similar domains
my $protein = ""; 

my $lek = "";

############################################
### 1. Read in the hmmsearch output file ###
############################################

print "\n*** reading in $hmmsearch_output ***\n\n";

@output = get_file_data("$hmmsearch_output");

print "*** Parsing $hmmsearch_output ***\n\n";

# If wrong format exit the program, ninth row from the end
if ($output [-9] !~ /^# Program:         hmmsearch.*/) {
	print "wrong file format. $hmmsearch_output has to be a hmmsearch output file\n\n";
	exit;
}

### 1.1 Trim input file (hmmsearch) to fit script. Erase lines starting with #. 
@output = grep !/^#.*/, @output;

# Read domain-specific coverage values file
# domname	cut-off
# EIN3		0.5874125874

open(FILENAME,$domspec_cuts);
my %cuts = map { chomp; split /\t/ } <FILENAME>;

### 1.2 Use first (prot name), third (domain name), eight (dom evalue), fourteenth (# overlapping) and fifteenth (# domain) column of data input(array element).	
foreach my $output (@output) {
	$output =~ /^(\S+)\s+\S+\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+.*/;
	$lek = (($6-$5)+1)/$3;
	$output = $1."\t".$2."\t".$4."\t".$lek;	
	}

# Stucture:
# ARATHwo_AT2G03430.1	Ank	1.8e-08		0.93
# ARATHwo_AT2G03430.1	Ank     7.8e-10		0.95
# ARATHwo_AT2G03430.1	Ank     1.2e-08		0.93
# ARATHwo_AT2G03430.1	Ank     1.2e-08		0.63

################################
### 2. Modify the input data ###
################################

### 2.1 Run through length_cut_off loop (throws out entries below coverage values)
my @cutarray;
foreach my $output (@output) {
	$output =~ /^(\S+)\t+(\S+)\t(\S+)\t(\S+)/;
	if ($4 > $cuts{$2}) {
	push (@cutarray, $output);
	}
}

# End up with (prot name) (domain name) (dom evalue) (overlapping)

# Structure:
# ARATHwo_AT2G03430.1	Ank	1.8e-08		0.93
# ARATHwo_AT2G03430.1	Ank     7.8e-10		0.95
# ARATHwo_AT2G03430.1	Ank     1.2e-08		0.93


### 2.2 Retain hit with the lowest evalue

# String for the previous protein name
my $old_prot = "";
# Coordinates for array filling
my $x = -1;
my $y = 1;
my $scoreeval;
my @newarray;
foreach my $cutarray (@cutarray) {
	$cutarray =~ /^(\S+)\t+(\S+)\t(\S+)\t(\S+)/;
	$akt_entry = $1.$2; 
	# If the protein contains more than one of the same domain
	if ($old_prot eq $akt_entry) {
		if ($3>$scoreeval) {
			$y++;
			$newarray[$x] = $1."\t".$2."\t".$3."\t".$y;
			}
		else {
			$y++;
			$newarray[$x] = $1."\t".$2."\t".$scoreeval."\t".$y;
			}	
		}
	# If the protein is new reset $j and push the entry in the next array
	else {
		$y = 1;
		$x++;
		$newarray[$x] = $1."\t".$2."\t".$3."\t".$y;
	}
	$old_prot = $akt_entry;
	$scoreeval = $3

}
# Final structure:
# ARATHwo_AT2G03430.1	Ank     7.8e-10		3

### 2.3 Then sort @output primarily on prot name($1) and then secondarily on eval($3) to fit the "similar"-procedure.
@newarray = sort { (split '\t', $a)[0] cmp (split '\t', $b)[0] || (split '\t', $a)[2] <=> (split '\t', $b)[2] } @newarray;

# Flags for merged domains
my $ARR_B_counter = 0;
my $bZIP_counter = 0;
my $ET_counter = 0;

# Flags for similar domain omission
my $WD40domain = 0;		# modified for FIE rule
my $FIEdomain = 0;		# --------'''''--------
my $MYBdomain = 0;
my $G2domain = 0;
my $YBdomain = 0;
my $YCdomain = 0;
my $PHDdomain = 0;
my $Alfindomain = 0;
my $Dr1domain = 0;
my $GATAdomain = 0;
my $Dofdomain = 0;
my $C1domain = 0;

foreach my $line (@newarray) {

# Reset the flags for domain omission
	if ($line !~ /^$protein.*$/ ) {
		$WD40domain = 0;
		$FIEdomain = 0;
		$MYBdomain = 0;
		$G2domain = 0;
		$YBdomain = 0;
		$YCdomain = 0;
		$PHDdomain = 0;
		$Alfindomain = 0;
		$Dr1domain = 0;
		$GATAdomain = 0;
		$Dofdomain = 0;
		$C1domain = 0;

	}
	
	# Get the Query entry
	$line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
	$protein = "$1";
	$extracted_domain = $2;
	$present = "$4"; # FIX! number for many Myb
	$whole_entry = $1.";".$2.";".$3.";".$4;


	# Ignore FIE if it is not better scored than WD40
	if ($extracted_domain eq "WD40" and $FIEdomain == 0) {$WD40domain = 1;}	 

	# Ignore similar domains part 1

	if ($extracted_domain eq "Myb_DNA-binding" and $G2domain == 0) {$MYBdomain = 1;}
	if ($extracted_domain eq "G2-like_Domain" and $MYBdomain == 0) {$G2domain = 1;}

	if ($extracted_domain eq "PHD" and $Alfindomain == 0 and $C1domain == 0) {$PHDdomain = 1;}
	if ($extracted_domain eq "Alfin-like" and $PHDdomain == 0 and $C1domain == 0) {$Alfindomain = 1;}
	if ($extracted_domain eq "C1_2" and $Alfindomain == 0 and $PHDdomain == 0) {$C1domain = 1;}

        if ($extracted_domain eq "GATA" and $Dofdomain == 0) {$GATAdomain = 1;}
	if ($extracted_domain eq "zf-Dof" and $GATAdomain == 0) {$Dofdomain = 1;}

	# Ignore similar domains part 2

	if ($extracted_domain eq "FIE_clipped_for_HMM" and $WD40domain == 1) {next;}
	if ($extracted_domain eq "Myb_DNA-binding" and $G2domain == 1) {next;}
	if ($extracted_domain eq "G2-like_Domain" and $MYBdomain == 1) {next;}

	if ($extracted_domain eq "PHD" and ($Alfindomain == 1 or $C1domain == 1)) {next;}
	if ($extracted_domain eq "Alfin-like" and ($PHDdomain == 1 or $C1domain ==1)) {next;}
	if ($extracted_domain eq "C1_2" and ($PHDdomain == 1 or $Alfindomain == 1)) {next;}

       	if ($extracted_domain eq "GATA" and $Dofdomain == 1) {next;}
	if ($extracted_domain eq "zf-Dof" and $GATAdomain == 1) {next;}
	
	# Save the entry
	push @results_of_extraction, $whole_entry;
	$entry_counter++;	
}
print "$entry_counter domain matches were found in $hmmsearch_output\n\n";

# Sort the entries:
my @sorted_results_of_extraction = @results_of_extraction;
print "*** Preparing hmmsearch results for classification ***\n\n";

# This array will be filled with the sorted proteins
my $array_of_arrays = [];

# String for the previous protein name
my $old_entry = "";

# Coordinates for array filling
my $i = -1;
my $j = 0;

@dec_table = get_file_data("$decision_table");

if ($dec_table [0] !~ /^[^;]+;[^;]+;[^;]/) {
	print "wrong file format. $decision_table has to be in the format family;domain;type\n\n";
	exit;
}

# Formate the entries
foreach my $domain_entry (@sorted_results_of_extraction) {

	$domain_entry =~ /^([^;]+);/;
	$akt_entry = $1; 
	# If the protein has more than one domain push it behind the old entry
	if ($old_entry eq $akt_entry) {
		$j++;
		$array_of_arrays->[$i][$j] = $domain_entry;
	}	

	# If the protein is new reset $j and push the entry in the next array
	else {
		$j = 0;
		$i++;
		$array_of_arrays->[$i][$j] = $domain_entry;
	}
	$old_entry = $akt_entry;

}
$old_entry = "";

# Remember the number of entries to use them later
my $number_of_proteins = $i+1;
print "$number_of_proteins proteins were found in $hmmsearch_output\n";
$number_of_proteins--;

#########################################
### 3. Get the classification entries ###
#########################################

# Make a list of all families in the classification rules file
my @liste_alle_familien = ('0_no_family_found');

my $ARR_B_in_list = 0;
my $bZIP_in_list = 0;
my $ET_in_list = 0;

my $array_of_classifications = [];
$i = -1;
$j = 0;

my @classifications_input = get_file_data("$decision_table");

# Sort the classifications
my @classifications = sort @classifications_input;


foreach my $classification (@classifications) {
	$classification =~ /^([^;]+);/;
	$akt_entry = $1;
	
	# If the family has more than one domain entry push it behind the old entry
	if ($old_entry eq $akt_entry) {
	    
		$j++;
		$array_of_classifications->[$i][$j] = $classification;
	}
			
	# If the family is new push it in the next array and reset $j
	else {
		$j = 0;
		$i++;
		$array_of_classifications->[$i][$j] = $classification;
		# Add family to list of all families in the classification rules
		push(@liste_alle_familien,$akt_entry);
		
		# Add hardcoded "OR" defined families to list of families if subfamily is present 
		# Do this just once for each family
		if (($akt_entry eq "GARP_ARR-B_Myb") or ($akt_entry eq "GARP_ARR-B_G2")) {
                	if ($ARR_B_in_list == 0) {
				push(@liste_alle_familien,"GARP_ARR-B");
				$ARR_B_in_list++;
			}}
		if (($akt_entry eq "bZIP1") or ($akt_entry eq "bZIP2") or ($akt_entry eq "bZIPAUREO") or ($akt_entry eq "bZIPCDD")) {
                	if ($bZIP_in_list == 0) {
				push(@liste_alle_familien,"bZIP");
				$bZIP_in_list++;
			}}
		if (($akt_entry eq "HRT") or ($akt_entry eq "GIY_YIG")) {
                	if ($ET_in_list == 0) {
				push(@liste_alle_familien,"ET");
				$ET_in_list++;
			}}
	}
	$old_entry = $akt_entry;
}	

my $number_of_classifications = $i+1;
print "$number_of_classifications different families were found in $decision_table\n\n"; 
$number_of_classifications--;

print "*** begin classification ***\n\n";

# Create the output files (output.1 and output.3)
my $outputfile = "$family_classifications";
my $subfamilyoutput = "$subfamily_classifications";

unless (open(FAMILY_CLASSIFICATIONS, ">$outputfile")) {
	print "Cannot open file \"$outputfile\" to write to!!\n\n";
	exit;
}

unless (open(SUBFAMILY_CLASSIFICATIONS, ">$subfamilyoutput")) {
	print "Cannot open file \"$subfamilyoutput\" to write to!!\n\n";
	exit;
}

print FAMILY_CLASSIFICATIONS "family classifications for $hmmsearch_output\n";
print SUBFAMILY_CLASSIFICATIONS "family classifications for $hmmsearch_output\n";

# Counter for classified families
my $classified_families = 0;
# Counter for unclassified families
my $unclassified_families = 0;
# Counter for the different amount of entries in classification possibilities
my @entry_counter = [];
# Counter for the possible classifications for the current protein
my $possibilities = 0;
# For later famliy statistics
my @family_list = [];
# This array will be filled with the result. Later written in output.1 and output.3
my @family_classifications_file = [];
my @subfamily_classifications_file = [];

# One time for every protein in array of arrays
for (my $ii = 0; $ii <= $number_of_proteins; $ii++) {
	my $jj = 0;
	$possibilities = 0;

    # Flag to mark that any family was found for the current protein
	my $member_found = 0;
	
	# In this string all the domains of one protein will be stored
	# The ; signs are important to mark a single entry (for example to differ AP2 and TF_AP2)
	my $domains = ";";

	# Get the domains in all elements of the current protein array:
	while ($array_of_arrays->[$ii][$jj]) {

		$array_of_arrays->[$ii][$jj] =~ /([^;]+);[^;]+;(\d+)$/;
		# For discrimination between MYB and MYB-related. If more than 1 MYB domains was found replaced "Myb_DNA-binding" with "two_or_more_Myb_DNA-binding"
		# Same for LIM
		if ($1 eq "Myb_DNA-binding" and $2 eq 2) {
			$domains .= "MYB-2R;Myb_DNA-binding;";
		}
		elsif ($1 eq "Myb_DNA-binding" and $2 eq 3) {
			$domains .= "MYB-3R;Myb_DNA-binding;";
		}
		elsif ($1 eq "Myb_DNA-binding" and $2 eq 4) {
			$domains .= "MYB-4R;Myb_DNA-binding;";
		}
		elsif ($1 eq "LIM" and $2 > 1) {
			$domains .= "two_or_more_LIM;LIM;";
		}
		else {
			$domains .= "$1;";
		}
		$jj++;
	}

#########################
### 4. Classification ###
#########################

# Get the name of the actual protein for the output after classification:
$array_of_arrays->[$ii][0] =~ /^([^;]+);/;
$akt_entry = $1;

# Variables for double classification decisions:

# Here the shoulds for the actual family will be stored
my $shoulds_act_family = ";";
# Here the shoulds of the previous CLASSIFIED family will be stored
my $shoulds_prev_family = "";

# This loop will be made for every family entry in output.1 and output.3
for ($i = 0;$i <= $number_of_classifications; $i++) {
	$j = 0;

	# Flag to mark if the protein is a possible member of the current family
	my $possible_member = 0;

	# Reset the domain entries for the current family
        $shoulds_act_family = ";";

        # This loop will be made for every domain classification entry for the current family
	while ($array_of_classifications->[$i][$j]) {
		
		# For "should" entries
		if ($array_of_classifications->[$i][$j] =~ /;should$/) {

			# Get the domain from the classification entry
			$array_of_classifications->[$i][$j] =~ /^[^;]+;([^;]+);[^;]+/;
			
			my $domain_classification_entry = $1;
			# Check if the protein might be a member of the current family. If yes check the next entry
			if ($domains =~ /;$domain_classification_entry;/) {
				$possible_member = 1;
				$j++;
				# All matched domains for the current family are stored here:
				$shoulds_act_family .= "$domain_classification_entry;";
				next;
			}
			$possible_member = 0;
			$j++;
			last;
		}

		# For "should not" entries
		if ($array_of_classifications->[$i][$j] =~ /;should not$/) {

			# Get the domain from the classification entry
			$array_of_classifications->[$i][$j] =~ /^[^;]+;([^;]+);[^;]+/;
			
			# Check if the protein has the domain. If yes go to the next family
			my $domain_check = $1;
			if ($domains =~ /;$domain_check;/) {
					$possible_member = 0;
					last;
			}
			# If not check the next entry
			else {
					$j++;
			    		next;
			}
		}
			
			# This happens when there is no "should" or "should not" entry at the right position
			print "error in classification entry $array_of_classifications->[$i][$j]\n"
	}
		       
	# If the check was successful print the found classification
	if ($possible_member) {

		$array_of_classifications->[$i][0] =~ /^([^;]+);/;
		my $currentFamily = $1;
		$possibilities++;
		# The double classified proteins solution:
		if ($possibilities > 1) {

			# Fill an array with the current domains in $domains
			my @domains_array = [];
			my $number_of_domains = ($domains =~ tr/;/;/)-1;
			my $tempdomains = $domains;
			for (my $domain_counter = 0;$domain_counter != $number_of_domains;$domain_counter++) {

				$tempdomains =~ /^;([^;]+);/;
				$domains_array[$domain_counter] = $1;
				# Remove the inserted domain from domains
				$tempdomains =~ s/^;[^;]+;/;/;
			}
			# Compare the two possible families (actual and previous family)
			my $array_counter = 0;
			while ($array_counter < $number_of_domains) {
				if ($shoulds_prev_family =~ /$domains_array[$array_counter]/) {
					if ($shoulds_act_family =~ /$domains_array[$array_counter]/) {

						# If actual domain is in previous and actual family go to next domain
						$array_counter++;
						next;

					}
					# Protein was right classified in the previous classification. The actual entry will be ignored

					$possibilities--;
					last;
				}

				# If actual domain is in the actual family but not in the previous. The actual classification is right
				if ($shoulds_act_family =~ /$domains_array[$array_counter]/) {

					pop @family_classifications_file;
					pop @subfamily_classifications_file;
					
					# Merge "GARBP_ARR-B_Myb" and "GARP_ARR-B_G2" to "GARP_ARR-B"
					if ( ($currentFamily eq "GARP_ARR-B_Myb") or ($currentFamily eq "GARP_ARR-B_G2")){

						push @family_classifications_file, "$akt_entry;GARP_ARR-B;$possibilities$domains\n";
						push @subfamily_classifications_file, "$akt_entry;GARP_ARR-B;-;$possibilities$domains\n";
						if ($ARR_B_counter == 0) {
							push @liste_alle_familien, 'GARP_ARR-B';
							$ARR_B_counter++;
						}
					}
					elsif ( ($currentFamily eq "bZIP1") or ($currentFamily eq "bZIP2") or ($currentFamily eq "bZIPAUREO") or ($currentFamily eq "bZIPCDD")) {

						push @family_classifications_file, "$akt_entry;bZIP;$possibilities$domains\n";
						push @subfamily_classifications_file, "$akt_entry;bZIP;-;$possibilities$domains\n";
						if ($bZIP_counter == 0) {
							push @liste_alle_familien, 'bZIP';
							$bZIP_counter++;
						}
					}
					elsif ( ($currentFamily eq "HRT") or ($currentFamily eq "GIY_YIG")) {
						push @subfamily_classifications_file, "$akt_entry;ET;-;$possibilities$domains\n";
						push @family_classifications_file, "$akt_entry;ET;$possibilities$domains\n";
						if ($ET_counter == 0) {
							push @liste_alle_familien, 'ET';
							$ET_counter++;
						}
					}
					else {

						push @family_classifications_file, "$akt_entry;$currentFamily;$possibilities$domains\n";
						push @subfamily_classifications_file, "$akt_entry;$currentFamily;-;$possibilities$domains\n";
					}	
					$possibilities--;
					last;

				}

					$array_counter++;
			}

		}
		# Normal enty: if the family is classified for one family write it in the array
		else {

			# Store the matched domains of the actual family for respectively later doubly classification decision 
			$shoulds_prev_family = $shoulds_act_family;
			# Store the entry in the output array
			# Merge "GARBP_ARR-B_Myb" and "GARP_ARR-B_G2" to "GARP_ARR-B"
			# Add GARP_ARR-B to @liste_alle_familien ONCE
				
				if ( ($currentFamily eq "GARP_ARR-B_Myb") or ($currentFamily eq "GARP_ARR-B_G2")){
					push @family_classifications_file, "$akt_entry;GARP_ARR-B;$possibilities$domains\n";
					push @subfamily_classifications_file, "$akt_entry;GARP_ARR-B;-;$possibilities$domains\n";
					if ($ARR_B_counter == 0) {
						push @liste_alle_familien, 'GARP_ARR-B';
						$ARR_B_counter++;
					}
				}
				elsif (($currentFamily eq "bZIP1") or ($currentFamily eq "bZIP2") or ($currentFamily eq "bZIPAUREO") or ($currentFamily eq "bZIPCDD")){
					push @family_classifications_file, "$akt_entry;bZIP;$possibilities$domains\n";
					push @subfamily_classifications_file, "$akt_entry;bZIP;-;$possibilities$domains\n";
					if ($bZIP_counter == 0) {
						push @liste_alle_familien, 'bZIP';
						$bZIP_counter++;
					}
				}
				elsif (($currentFamily eq "HRT") or ($currentFamily eq "GIY_YIG")){
					push @family_classifications_file, "$akt_entry;ET;$possibilities$domains\n";
					push @subfamily_classifications_file, "$akt_entry;ET;-;$possibilities$domains\n";
					if ($ET_counter == 0) {
						push @liste_alle_familien, 'ET';
						$ET_counter++;
					}
				}
				else {

				push @family_classifications_file, "$akt_entry;$currentFamily;$possibilities$domains\n";
				push @subfamily_classifications_file, "$akt_entry;$currentFamily;-;$possibilities$domains\n";
				}

                	$classified_families++;
			}

		# Flag to mark that any family was found for the current protein
		$member_found = 1;
		}

	}
	
	if ($possibilities == 0) {

		# If no family was found for the current protein give this out
		push @family_classifications_file, "$akt_entry;0_no_family_found;0$domains\n";
		push @subfamily_classifications_file, "$akt_entry;0_no_family_found;-;0$domains\n";
		$unclassified_families++;
	}
	
	# Increment the number of possible members in the array at the right position
	else {

		my $temp = $entry_counter[$possibilities];
		$temp++;
		$entry_counter[$possibilities] = $temp;
	}
	
}
# Define for which TAP families also subfamilies may be present:

foreach my $subfamily_classifications_file (@subfamily_classifications_file) {
	
	# If "C2H2" was assigned to a sequence, write "C2H2;C2H2" to output.3 as TAP family and subfamily, respectively. 
	if ($subfamily_classifications_file =~ /C2H2/ ) {
		$subfamily_classifications_file =~ s/C2H2;-/C2H2;C2H2/
	}
	# If "C2H2_IDD" was assigned to a sequence, write "C2H2;C2H2_IDD" to output.3 as TAP family and subfamily, respectively. 
	if ($subfamily_classifications_file =~ /C2H2_IDD/ ) {
		$subfamily_classifications_file =~ s/C2H2_IDD;-/C2H2;C2H2_IDD/
	}
	if ($subfamily_classifications_file =~ /bHLH/ ) {
		$subfamily_classifications_file =~ s/bHLH;-/bHLH;bHLH/
	}
	if ($subfamily_classifications_file =~ /bHLH_TCP/ ) {
		$subfamily_classifications_file =~ s/bHLH_TCP;-/bHLH;bHLH_TCP/
	}
	if ($subfamily_classifications_file =~ /NF-YA/ ) {
		$subfamily_classifications_file =~ s/NF-YA;-/NFY;NF-YA/
	}
	if ($subfamily_classifications_file =~ /NF-YB/ ) {
		$subfamily_classifications_file =~ s/NF-YB;-/NFY;NF-YB/
	}
	if ($subfamily_classifications_file =~ /NF-YC/ ) {
		$subfamily_classifications_file =~ s/NF-YC;-/NFY;NF-YC/
	}
	if ($subfamily_classifications_file =~ /C1HDZ/ ) {
		$subfamily_classifications_file =~ s/C1HDZ;-/HDZ;C1HDZ/
	}
	if ($subfamily_classifications_file =~ /C2HDZ/ ) {
		$subfamily_classifications_file =~ s/C2HDZ;-/HDZ;C2HDZ/
	}
	if ($subfamily_classifications_file =~ /C3HDZ/ ) {
		$subfamily_classifications_file =~ s/C3HDZ;-/HDZ;C3HDZ/
	}
	if ($subfamily_classifications_file =~ /C4HDZ/ ) {
		$subfamily_classifications_file =~ s/C4HDZ;-/HDZ;C4HDZ/
	}
	if ($subfamily_classifications_file =~ /MYST/ ) {
		$subfamily_classifications_file =~ s/MYST;-/HAT;MYST/
	}
	if ($subfamily_classifications_file =~ /CBP/ ) {
		$subfamily_classifications_file =~ s/CBP;-/HAT;CBP/
	}
	if ($subfamily_classifications_file =~ /TAFII250/ ) {
		$subfamily_classifications_file =~ s/TAFII250;-/HAT;TAFII250/
	}
	if ($subfamily_classifications_file =~ /GNAT/ ) {
		$subfamily_classifications_file =~ s/GNAT;-/HAT;GNAT/
	}
	if ($subfamily_classifications_file =~ /LOB1/ ) {
		$subfamily_classifications_file =~ s/LOB1;-/LBD;LOB1/
	}
	if ($subfamily_classifications_file =~ /LOB2/ ) {
		$subfamily_classifications_file =~ s/LOB2;-/LBD;LOB2/
	}
	if ($subfamily_classifications_file =~ /MYB-related/ ) {
		$subfamily_classifications_file =~ s/MYB-related;-/MYB;MYB-related/
	}
	if ($subfamily_classifications_file =~ /SWI\/SNF_SWI3/ ) {
		$subfamily_classifications_file =~ s/SWI\/SNF_SWI3;-/MYB-related;SWI\/SNF_SWI3/
	}
	if ($subfamily_classifications_file =~ /MYB-2R/ ) {
		$subfamily_classifications_file =~ s/MYB-2R;-/MYB;MYB-2R/
	}
	if ($subfamily_classifications_file =~ /MYB-3R/ ) {
		$subfamily_classifications_file =~ s/MYB-3R;-/MYB;MYB-3R/
	}
	if ($subfamily_classifications_file =~ /MYB-4R/ ) {
		$subfamily_classifications_file =~ s/MYB-4R;-/MYB;MYB-4R/
	}
	if ($subfamily_classifications_file =~ /RKD/ ) {
		$subfamily_classifications_file =~ s/RKD;-/RWP-RK;RKD/
	}
	if ($subfamily_classifications_file =~ /NLP/ ) {
		$subfamily_classifications_file =~ s/NLP;-/RWP-RK;NLP/
	}
	if ($subfamily_classifications_file =~ /AP2/ ) {
		$subfamily_classifications_file =~ s/AP2;-/AP2;AP2/
	}
	if ($subfamily_classifications_file =~ /CRF/ ) {
		$subfamily_classifications_file =~ s/CRF;-/AP2;CRF/
	}
}

### 4.1 Filter for gene models

# Filter for gene models. All alternative splice variants will be removed from Arabidopsis and rice entries
shift @family_classifications_file;
shift @subfamily_classifications_file;
# Sort the array for filter and a sorted output
@family_classifications_file = sort @family_classifications_file;
@subfamily_classifications_file = sort @subfamily_classifications_file;
if ($gene_model_filter and $gene_model_filter eq "filter") {
	print "*** filtering gene models: removing splice variants ***\n";
	
	# temporary array for the filtered entries
	my @filtered_entries = [];
	my $models_removed_counter = 0;

	my $prev_entry = "";
	foreach my $entry_line (@family_classifications_file) {
		if ($entry_line !~ /^([^.]+).\d+/) {
			print "Bad format for filtering splice variants. Use \"filter\" only for TAIR (Aradopsis) or TIGR (Rice) files.\n\n";
			exit;
		}
		if ($1 eq $prev_entry) {
			$models_removed_counter++;
			next;
		}
		$prev_entry = $1;
		push @filtered_entries, $entry_line;
	}
	print " -> $models_removed_counter gene models were removed\n\n";
	@family_classifications_file = @filtered_entries;
	shift @family_classifications_file;
}

# Print the results in the array to output.1 and output.3 and push the names of the TAP families and Subfamilies into $family_list to create output.2
foreach my $fcf_line (@family_classifications_file) {
	$fcf_line =~ /^[^;]+;([^;]+)/;
	#push @family_list, "$1";
	print FAMILY_CLASSIFICATIONS "$fcf_line";
}
close (FAMILY_CLASSIFICATIONS);

foreach my $fcf_line (@subfamily_classifications_file) {
	$fcf_line =~ /^[^;]+;([^;]+);([^;]+)/;
	if ($2 eq "-") {
		push @family_list, "$1";
	}
	else {
		push @family_list, "$2";
	}
	print SUBFAMILY_CLASSIFICATIONS "$fcf_line";
} 

close (SUBFAMILY_CLASSIFICATIONS);

print "*** calculating the family statistics and write it in $family_statistics ***\n\n";

##################################
### 5. Create the output files ###
##################################

my $statistics_outputfile = "$family_statistics";

unless (open(FAMILY_STATISTICS, ">$statistics_outputfile")) {
	print "Cannot open file \"$statistics_outputfile\" to write to!!\n\n";
	exit;
}

# Count the family entries
my @output_family_statistics = ();
my @gefundene_familien = ();
my $family_counter = 1;

shift @family_list;
@family_list = sort @family_list;


my $old_family = "";
push @family_list, 'BAD FIX'; # Makes to loop go through every fam and stops at non fam(BAD FIX).

foreach my $line (@family_list) {
	if ($line eq $old_family) {
		$family_counter++;
	}
	elsif ($old_family ne "") {
		push (@output_family_statistics,"$old_family;$family_counter\n");
		# Add all found families to a list of found families
		push (@gefundene_familien,"$old_family");
	        $family_counter=1;
	}
	$old_family = $line;
}

my %hash = ();

# Put all families from the classifictaion ruled and all found families in a hash
foreach my $element (@gefundene_familien,@liste_alle_familien) {$hash{$element}++;}

# Remove merged families
delete $hash{'GARP_ARR-B_Myb'};
delete $hash{'GARP_ARR-B_G2'};
delete $hash{'bZIP1'};
delete $hash{'bZIP2'};
delete $hash{'bZIPAUREO'};
delete $hash{'bZIPCDD'};
delete $hash{'HRT'};
delete $hash{'GIY_YIG'};

# Add all not found families from the classification rules to @output_family_statistics 
# With zero as number of families found 

foreach my $element (keys %hash) {
	if (($hash{$element} == 1) and ( $element eq "0_no_family_found")) {
	push (@output_family_statistics,"$element;$unclassified_families\n");
	}
	if (($hash{$element} == 1) and ($element ne "0_no_family_found")) {
	push (@output_family_statistics,"$element;0\n");
	}
}

# Sort @output_family_statistics caseinsensitive-alphabetically 
my @sortierte_statistik = sort {lc $a cmp lc $b} @output_family_statistics;

# Print headline to @sortierte_statistik
unshift (@sortierte_statistik,"family statistics for $hmmsearch_output\n");

print FAMILY_STATISTICS @sortierte_statistik;

# Print FAMILY_STATISTICS "$old_family;$family_counter\n";

close (FAMILY_STATISTICS);

#################################################
### 6. Give out some statistical informations ###
#################################################

$entry_counter [0] = 0;
my $sum = 0;
foreach my $entry (@entry_counter) {
	$sum += $entry;
}
print "$classified_families classifications were found for $sum proteins.\n";
print "This classifications are divided in:\n";
my $count = 0;
foreach my $element (@entry_counter) {
	if ($count != 0) {
		print "$element proteins were classified for $count";
		if ($count == 1) {print " family\n";}
		else {print " different families\n";}
	}
	$count++;
}
print "\n$unclassified_families proteins could not be classified\n\n";

print "*** The results were written in $family_classifications and $subfamily_classifications ***\n";
print "*** done ***\n\n";

exit;

sub get_file_data {
	
	my ($filename) = @_;

	use strict;
	use warnings;

	my @filedata = ();

	unless( open(GET_FILE_DATA, $filename)) {
		print STDERR "Cannot open file \"$filename\"n\n";
		exit;
	}

	@filedata = <GET_FILE_DATA>;

	close GET_FILE_DATA;

	return @filedata;
}

