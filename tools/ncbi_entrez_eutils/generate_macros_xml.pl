#!/usr/bin/env perl

#Usage: perl generate_macros_xml.pl > macros.xml

#Note, this script uses einfo.py to get database info.  It also uses manually compiled data stored at the bottom of this script that is based on: https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
#The data in the table on that page was manipulated to replace nulls with 'none', remove duplicates, and add missing formats based on correspondence with MLN.

##
## use einfo to retrieve all the valid databases
##

print STDERR "Retrieving database list\n";

my $dbxml = `python einfo.py --user_email "planemo@galaxyproject.org" --admin_email "planemo@galaxyproject.org;test@bx.psu.edu"`;

my(@dblist);
my $dbs     = {};
my $dbfroms = {};
my $dbnames = {};
foreach(split(/\n/,$dbxml))
  {
    if(/<DbName>(.+)<\/DbName>/)
      {
        my $db = $1;
        push(@dblist,$db);
        $dbs->{$db}     = 0;
        $dbfroms->{$db} = 0;
        $dbnames->{$db} = $_;
      }
  }

##
## Use einfo to retrieve all the valid links for each database (Note: some databases are not linked)
##

my $h = {};
foreach my $db (sort {$dbnames->{$a} cmp $dbnames->{$b}} @dblist)
  {
    sleep(2);

    print STDERR "Retrieving info for $db\n";

    my $response = `python einfo.py --db $db --user_email "planemo\@galaxyproject.org" --admin_email "planemo\@galaxyproject.org;test\@bx.psu.edu"`;

    my $dolinks = 0;
    my $link    = "";
    my $name    = "";

    foreach(split(/\n/,$response))
      {
        if(/<LinkList>/)
          {
            $dolinks = 1;
            #Save whether there exist links from this database
            $dbfroms->{$db} = 1;
          }
        elsif(!$dolinks)
          {
            if(/<MenuName>(.+)<\/MenuName>/)
              {$dbnames->{$db} = "$1 ($db)"}
          }
        elsif($dolinks)
          {
            if(/<Name>(.+)<\/Name>/)
              {$link=$1}
            elsif(/<Menu>(.*)<\/Menu>/)
              {$name=$1}
            elsif(/<DbTo>(.+)<\/DbTo>/)
              {
                $dbto=$1;
                push(@{$h->{$db}->{$dbto}},[$link,$name]);
                $link="";
                $name="";
              }
          }
      }
  }

my @sorted_dblist = sort {$dbnames->{$a} cmp $dbnames->{$b}} @dblist;

##
## Generate XML to govern the valid databases to use with efetch
##

my $efetch_dbhash = {}; #->{efetch-compatible-db}->{rettype-retmode-galaxy_format} = format_name (galaxy_format)
while(<DATA>)
  {
    chomp;
    my($db,$galaxy_format,$retmode,$rettype,$format_name) = split(/\t/,$_);
    $efetch_dbhash->{$db}->{"$rettype-$retmode-$galaxy_format"} =
      "$format_name ($galaxy_format)";
  }

#EFetch database select list

print << 'EOXML';
  <xml name="dbselect_efetch" token_name="db_select" token_label="NCBI Database to Query">
    <param name="@NAME@" type="select" label="@LABEL@">
EOXML

foreach my $db (grep {exists($dbs->{$_})}
                sort {$dbnames->{$a} cmp $dbnames->{$b}}
                keys(%$efetch_dbhash))
  {
    my $selected = '';
    if($db eq 'pubmed')
      {$selected = ' selected="True"'}
    print << "    EOXML";
      <option value="$db"$selected>$dbnames->{$db}</option>
    EOXML
  }

print << 'EOXML';
    </param>
  </xml>
EOXML

#EFetch output formats

print << 'EOXML';
  <xml name="efetchdb">
    <conditional name="db">
      <expand macro="dbselect_efetch" />
EOXML

foreach my $db (grep {exists($dbs->{$_})}
                sort {$dbnames->{$a} cmp $dbnames->{$b}}
                keys(%$efetch_dbhash))
  {
    print << "    EOXML";
      <when value="$db">
        <param name="output_format" type="select" label="Output Format">
    EOXML

    foreach my $eutils_format (sort {$efetch_dbhash->{$db}->{$a} cmp
                                       $efetch_dbhash->{$db}->{$b}}
                               keys(%{$efetch_dbhash->{$db}}))
      {
        print << "        EOXML";
          <option value="$eutils_format">$efetch_dbhash->{$db}->{$eutils_format}</option>
        EOXML
      }

    print << "    EOXML";
        </param>
      </when>
    EOXML
  }

print << 'EOXML';
    </conditional>
  </xml>
EOXML

##
## Create a select list for the databases linked *from*
##

print << 'EOXML';
  <xml name="dbselect" token_name="db_select" token_label="NCBI Database to Query">
    <param name="@NAME@" type="select" label="@LABEL@">
EOXML

foreach my $from (@sorted_dblist)
  {
    print << "    EOXML";
      <option value="$from">$dbnames->{$from}</option>
    EOXML
  }

print << 'EOXML';
    </param>
  </xml>
EOXML

##
## Create a select list for the databases linked *to*
##

print << 'EOXML';
  <xml name="dbselect_linked" token_name="db_select_linked" token_label="NCBI Database to Use">
    <param name="@NAME@" type="select" label="@LABEL@">
EOXML

foreach my $from (grep {$dbfroms->{$_}} @sorted_dblist)
  {
    print << "    EOXML";
      <option value="$from">$dbnames->{$from}</option>
    EOXML
  }

print << 'EOXML';
    </param>
  </xml>
EOXML

##
## Create empty entries for commands that take no *to* database or link
##

print << 'EOXML';
  <xml name="none_link_macro">
            <conditional name="db_to">
              <param name="db_select_to" type="select" label="To NCBI Database (n/a)">
                <option value="n/a">Not applicable</option>
              </param>
              <when value="n/a">
                <param name="linkname" type="select" label="Link Name (n/a)">
                  <option value="n/a">Not applicable</option>
                </param>
              </when>
            </conditional>
  </xml>
  <xml name="db_link_macro">
        <conditional name="db_from_link">
          <expand macro="dbselect_linked" name="db_select_from_link" label="From NCBI Database" />
EOXML

foreach(grep {$dbfroms->{$_}} @sorted_dblist)
  {
    print << "    EOXML";
          <when value="$_">
            <expand macro="none_link_macro" name="db_select_none" label="To NCBI Database" />
          </when>
    EOXML
  }

print << 'EOXML';
        </conditional>
  </xml>
EOXML

##
## This is the master macro for the command selection
##

print << 'EOXML';
  <xml name="linkmacro">
    <conditional name="cmd">
      <param name="cmd_select" type="select" label="Link Method" help="Fetch UIDs from the 'To' Database that are linked to supplied UIDs in the 'From' database">
        <option value="neighbor" selected="true">Neighbor (neighbor)</option>
        <option value="neighbor_history">Neighbor, save result in history server (neighbor_history)</option>
        <option value="neighbor_score">Neighbor Score (neighbor_score)</option>
        <option value="acheck">Show available links to any database (acheck)</option>
        <option value="ncheck">Show available links within the same database (ncheck)</option>
        <option value="lcheck">Show available links to external sources (LinkOuts) (lcheck)</option>
        <option value="llinks">Show available URLs and attributes for non-library LinkOut providers (llinks)</option>
        <option value="llinkslib">Show available URLs and attributes for all LinkOut Providers (llinkslib)</option>
        <option value="prlinks">Show available primary LinkOut Provider Links (prlinks)</option>
      </param>
      <when value="neighbor">
        <expand macro="db_db_link_macro" name="link_select" label="Link name" />
        <param name="output_format" type="select" label="Output Format">
          <option value="xml">ID File (xml)</option>
          <option value="json">ID File (json)</option>
          <option value="text" selected="true">ID File (tabular)</option>
        </param>
      </when>
      <when value="neighbor_history">
        <expand macro="db_db_link_macro" name="link_select" label="Link name" />
        <param name="output_format" type="select" label="Output Format">
          <option value="json">History File (json)</option>
          <option value="xml" selected="true">History File (xml)</option>
        </param>
      </when>
      <when value="neighbor_score">
        <expand macro="db_db_link_macro" name="link_select" label="Link name" />
        <param name="output_format" type="select" label="Output Format">
          <option value="xml">ID File (xml)</option>
          <option value="json">ID File (json)</option>
          <option value="text" selected="true">ID File (tabular)</option>
        </param>
      </when>
      <when value="acheck">
        <expand macro="db_link_macro" name="db_select_from_link" label="From NCBI Database" />
        <param name="output_format" type="select" label="Output Format">
          <option value="xml" selected="True">Link Description File (xml)</option>
          <option value="json">Link Description File (json)</option>
        </param>
      </when>
      <when value="ncheck">
        <expand macro="db_link_macro" name="db_select_from_link" label="From NCBI Database" />
        <param name="output_format" type="select" label="Output Format">
          <option value="xml" selected="True">Link Description File (xml)</option>
          <option value="json">Link Description File (json)</option>
        </param>
      </when>
      <when value="lcheck">
        <expand macro="db_link_macro" name="db_select_from_link" label="From NCBI Database" />
        <param name="output_format" type="select" label="Output Format">
          <option value="xml" selected="True">Link Description File (xml)</option>
          <option value="json">Link Description File (json)</option>
        </param>
      </when>
      <when value="llinks">
        <expand macro="db_link_macro" name="db_select_from_link" label="From NCBI Database" />
        <param name="output_format" type="select" label="Output Format">
          <option value="xml" selected="True">Link Description File (xml)</option>
          <option value="json">Link Description File (json)</option>
        </param>
      </when>
      <when value="llinkslib">
        <expand macro="db_link_macro" name="db_select_from_link" label="From NCBI Database" />
        <param name="output_format" type="select" label="Output Format">
          <option value="xml" selected="true">Link Description File (xml)</option>
          <option value="json">Link Description File (json)</option>
        </param>
      </when>
      <when value="prlinks">
        <expand macro="db_link_macro" name="db_select_from_link" label="From NCBI Database" />
        <param name="output_format" type="select" label="Output Format">
          <option value="xml" selected="true">Link Description File (xml)</option>
          <option value="json">Link Description File (json)</option>
        </param>
      </when>
    </conditional>
  </xml>
EOXML

##
## Create selections for valid links for command types neighbor, neighbor_history, and neighbor_score
##

print << 'EOXML';
  <xml name="db_db_link_macro">
        <conditional name="db_from_link">
          <expand macro="dbselect_linked" name="db_select_from_link" label="From NCBI Database" />
EOXML

foreach my $from (grep {$dbfroms->{$_}} @sorted_dblist)
  {
    print STDERR ("Creating Links From: $from\n");

    print << "    EOXML";
          <when value="$from">
            <conditional name="db_to">
              <param name="db_select_to" type="select" label="To NCBI Database">
    EOXML

    my @dbtos = (grep {exists($h->{$from}) && exists($h->{$from}->{$_})}
                 @sorted_dblist);
    foreach(@dbtos)
      {
        print << "        EOXML";
                <option value="$_">$dbnames->{$_}</option>
        EOXML
      }
    if(scalar(@dbtos) == 0)
      {
        #Provide an option for a self-link: from->from
        print << "        EOXML";
                <option value="$from">$dbnames->{$from}</option>
        EOXML
      }

    print << '    EOXML';
              </param>
    EOXML

    if(exists($h->{$from}))
      {
        #There do exist links to invalid(/outdated/non-existant) databases that
        #would result in an error if they are selected, so we use the original
        #@dblist instead of the keys present in the sub hash of $h->{$from}, and
        #then check for existence in the sub-hash
        foreach my $to (grep {exists($h->{$from}->{$_})} @sorted_dblist)
          {
            print STDERR ("\tTo: $to Links: ",
                          join(',',map {$_->[0]} @{$h->{$from}->{$to}}),
                          "\n");

            print << "            EOXML";
              <when value="$to">
                <param name="linkname" type="select" label="Link Name">
                  <option value="None">All Links</option>
            EOXML

            foreach(sort {"$a->[1] ($a->[0])" cmp "$b->[1] ($b->[0])"}
                    @{$h->{$from}->{$to}})
              {
                print << "                EOXML";
                  <option value="$_->[0]">$_->[1] ($_->[0])</option>
                EOXML
              }

            print << "            EOXML";
                </param>
              </when>
            EOXML

          }
      }
    else
      {
        ##
        ## Add-on selections for self-links for command types neighbor,
        ## neighbor_history, and neighbor_score
        ## Note, I'm not sure this would yield a valid result from elink
        ##

        #This shows $from, but this is the 'when' for db_to conditional
        print << "        EOXML";
              <when value="$from">
                <param name="linkname" type="select" label="Link Name">
                  <option value="none">All Links</option>
                </param>
              </when>
        EOXML
      }

    print << '    EOXML';
            </conditional>
          </when>
    EOXML
  }

##
## Add-on selections for self-links for command types neighbor,
## neighbor_history, and neighbor_score
## Note, I'm not sure this would yield a valid result from elink
##

foreach my $from (grep {!exists($h->{$_})} @sorted_dblist)
  {
    print << "EOXML";
          <when value=\"$from\">
            <conditional name=\"db_to\">
              <param name=\"db_select_to\" type=\"select\" label=\"To NCBI Database\">
                <option value=\"none\">Not applicable</option>
              </param>
              <when value=\"none\">
                <param name=\"linkname\" type=\"select\" label=\"Link Name\">
                  <option value=\"none\">Not applicable</option>
                </param>
              </when>
            </conditional>
          </when>
EOXML
  }

##
## This is the corresponding code for using the selections to add the respective command line options
##

print << 'EOXML';
    </conditional>
  </xml>
EOXML

print << 'EOXML';
  <token name="@LINK_TOKEN@">
    <![CDATA[
#if $cmd.db_from_link.db_to.db_select_to == 'n/a':
    none
#else:
    $cmd.db_from_link.db_to.db_select_to
#end if

$cmd.db_from_link.db_select_from_link

$cmd.cmd_select

#if $cmd.output_format == 'json':
    --retmode json
#elif $cmd.output_format == 'text':
    --retmode uilist
#else:
    --retmode xml
#end if

#if $cmd.db_from_link.db_to.linkname != 'None' and $cmd.cmd_select in ('neighbor', 'neighbor_history', 'neighbor_score'):
    --linkname $cmd.db_from_link.db_to.linkname
#end if
    ]]>
  </token>
EOXML

sub startXML
  {
    print << '    EOXML';
<?xml version="1.0"?>
<macros>
  <token name="@PROFILE@">18.01</token>
  <token name="@WRAPPER_VERSION@">1.70</token>
  <token name="@EMAIL_ARGUMENTS@">
--user_email "$__user_email__"
#set admin_emails = ';'.join(str($__admin_users__).split(','))
--admin_email "$admin_emails"
  </token>
  <!--  TODO: citation -->
  <token name="@REFERENCES@"><![CDATA[
  ]]></token>
  <token name="@DISCLAIMER@"><![CDATA[
Usage Guidelines and Requirements
=================================

Frequency, Timing, and Registration of E-utility URL Requests
-------------------------------------------------------------

In order not to overload the E-utility servers, NCBI recommends that users
limit large jobs to either weekends or between 9:00 PM and 5:00 AM Eastern time
during weekdays. Failure to comply with this policy may result in an IP address
being blocked from accessing NCBI.

Minimizing the Number of Requests
---------------------------------

If a task requires searching for and/or downloading a large number of
records, it is much more efficient to use the Entrez History to upload
and/or retrieve these records in batches rather than using separate
requests for each record. Please refer to Application 3 in Chapter 3
for an example. Many thousands of IDs can be uploaded using a single
EPost request, and several hundred records can be downloaded using one
EFetch request.


Disclaimer and Copyright Issues
-------------------------------

In accordance with requirements of NCBI's E-Utilities, we must provide
the following disclaimer:

Please note that abstracts in PubMed may incorporate material that may
be protected by U.S. and foreign copyright laws. All persons
reproducing, redistributing, or making commercial use of this
information are expected to adhere to the terms and conditions asserted
by the copyright holder. Transmission or reproduction of protected
items beyond that allowed by fair use (PDF) as defined in the copyright
laws requires the written permission of the copyright owners. NLM
provides no legal advice concerning distribution of copyrighted
materials. Please consult your legal counsel. If you wish to do a large
data mining project on PubMed data, you can enter into a licensing
agreement and lease the data for free from NLM. For more information on
this please see `https://www.nlm.nih.gov/databases/download/data_distrib_main.html <https://www.nlm.nih.gov/databases/download/data_distrib_main.html>`__

The `full disclaimer <https://www.ncbi.nlm.nih.gov/home/about/policies/>`__ is available on
their website

Liability
~~~~~~~~~

For documents and software available from this server, the
U.S. Government does not warrant or assume any legal liability or
responsibility for the accuracy, completeness, or usefulness of any
information, apparatus, product, or process disclosed.

Endorsement
~~~~~~~~~~~

NCBI does not endorse or recommend any commercial
products, processes, or services. The views and opinions of authors
expressed on NCBI's Web sites do not necessarily state or reflect those
of the U.S. Government, and they may not be used for advertising or
product endorsement purposes.

External Links
~~~~~~~~~~~~~~

Some NCBI Web pages may provide links to other Internet
sites for the convenience of users. NCBI is not responsible for the
availability or content of these external sites, nor does NCBI endorse,
warrant, or guarantee the products, services, or information described
or offered at these other Internet sites. Users cannot assume that the
external sites will abide by the same Privacy Policy to which NCBI
adheres. It is the responsibility of the user to examine the copyright
and licensing restrictions of linked pages and to secure all necessary
permissions.
        ]]></token>
  <token name="@LIST_OR_HIST@">
#if $query_source.qss == "history_json":
    --history_file $query_source.history_file
#else if $query_source.qss == "history_xml":
    --history_xml $query_source.history_xml
#else if $query_source.qss == "id_file":
    --id_list $query_source.id_file
#else if $query_source.qss == "id_list":
    --id $query_source.id_list
#else if $query_source.qss == "id_xml":
    --id_xml $query_source.id_xml
#else if $query_source.qss == "id_json":
    --id_json $query_source.id_json
#end if
  </token>
  <xml name="list_or_hist">
    <conditional name="query_source">
      <param name="qss" type="select" label="Enter Query IDs by..." help="Files output by ELink or ESearch are acceptable.  Query IDs in an ELink result are ignored.">
        <option value="history_json">History File (JSON)</option>
        <option value="history_xml">History File (XML)</option>
        <option value="id_file" selected="True">ID file (Tabular)</option>
        <option value="id_xml">ID File (XML)</option>
        <option value="id_json">ID File (JSON)</option>
        <option value="id_list">Paste IDs</option>
      </param>
      <when value="history_json">
        <param label="History File (JSON)" name="history_file" type="data" format="json" help="A JSON file containing the WebEnv ID and Query Key referencing the search on the NCBI history server"/>
      </when>
      <when value="history_xml">
        <param label="History File (XML)" name="history_xml" type="data" format="xml" help="An XML file containing the WebEnv ID and Query Key referencing the search on the NCBI history server"/>
      </when>
      <when value="id_file">
        <param label="ID File (Text)" name="id_file" type="data" format="txt,tabular" help="A Text file containing one ID per line"/>
      </when>
      <when value="id_xml">
        <param label="ID File (XML)" name="id_xml" type="data" format="xml" help="ESearch or ELink Result XML file"/>
      </when>
      <when value="id_json">
        <param label="ID File (JSON)" name="id_json" type="data" format="json" help="ESearch or ELink Result JSON file"/>
      </when>
      <when value="id_list">
        <param label="Paste ID List" name="id_list" type="text" area="true" help="Newline/Comma separated list of IDs"/>
      </when>
    </conditional>
  </xml>
  <xml name="citations">
    <citations>
      <citation type="bibtex">@Book{ncbiEutils,
          author = {Eric Sayers},
          title = {Entrez Programming Utilities Help},
          year = {2010},
          publisher = {National Center for Biotechnology Information, Bethesda, Maryland},
          note = {https://www.ncbi.nlm.nih.gov/books/NBK25500/}
      }</citation>
    </citations>
  </xml>
  <xml name="requirements">
    <requirements>
      <requirement type="package" version="1.70">biopython</requirement>
    </requirements>
  </xml>
  <token name="@EFETCH_FORMAT_TOKEN@">
    <![CDATA[

    ## This token must go at the end of the efetch command

    #set rettype, retmode, format = str($db.output_format).split('-')

    #if retmode != "none":
      --retmode $retmode
    #end if
    ## Otherwise, defaults to a None/empty which implies 'default' to NCBI

    #if rettype != "none":
        --rettype $rettype
    #end if

    --galaxy_format $format

    ]]>
  </token>
    EOXML
  }

sub endXML
  {
    print << '    EOXML';
</macros>
    EOXML
  }

BEGIN {startXML()}
END   {endXML()}


##
## Output formats for efetch mapped to galaxy formats
##

#Based on:
#https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly

#Note: While json works for esearch and elink, the only database that supports
#json (according to an NLM support ticket I have about this) is snp

#The output_format param value for these will be "rettype-retmode-format"

#db	galaxy	retmode	rettype	format_name
__DATA__
bioproject	tabular	text	uilist	List of UIDs
bioproject	xml	xml	docsum	Document summary
bioproject	xml	xml	uilist	List of UIDs
bioproject	xml	xml	xml	Full record
biosample	tabular	text	uilist	List of UIDs
biosample	txt	text	full	Full record
biosample	xml	xml	docsum	Document summary
biosample	xml	xml	full	Full record
biosample	xml	xml	uilist	List of UIDs
biosystems	tabular	text	uilist	List of UIDs
biosystems	xml	xml	docsum	Document summary
biosystems	xml	xml	uilist	List of UIDs
biosystems	xml	xml	xml	Full record
clinvar	tabular	text	uilist	List of UIDs
clinvar	xml	xml	clinvarset	ClinVar Set
clinvar	xml	xml	docsum	Document summary
clinvar	xml	xml	uilist	List of UIDs
clinvar	xml	none	none	Full
gds	tabular	text	uilist	List of UIDs
gds	txt	text	summary	Summary
gds	xml	xml	docsum	Document summary
gds	xml	xml	uilist	List of UIDs
gds	xml	none	none	Full
gene	txt	text	gene_table	Gene table
gene	tabular	text	uilist	List of UIDs
gene	txt	asn.1	none	text ASN.1
gene	xml	xml	docsum	Document summary
gene	xml	xml	none	Full
gene	xml	xml	uilist	List of UIDs
gtr	tabular	text	uilist	List of UIDs
gtr	xml	xml	docsum	Document summary
gtr	xml	xml	gtracc	GTR Test Report
gtr	xml	xml	uilist	List of UIDs
gtr	xml	none	none	Full
homologene	fasta	text	fasta	FASTA
homologene	tabular	text	alignmentscores	Alignment scores
homologene	tabular	text	uilist	List of UIDs
homologene	txt	asn.1	none	text ASN.1
homologene	txt	text	homologene	HomoloGene
homologene	xml	xml	docsum	Document summary
homologene	xml	xml	none	Full
homologene	xml	xml	uilist	List of UIDs
mesh	tabular	text	uilist	List of UIDs
mesh	txt	text	full	Full record
mesh	xml	xml	docsum	Document summary
mesh	xml	xml	uilist	List of UIDs
nlmcatalog	tabular	text	uilist	List of UIDs
nlmcatalog	txt	text	none	Full record
nlmcatalog	xml	xml	docsum	Document summary
nlmcatalog	xml	xml	none	Full
nlmcatalog	xml	xml	uilist	List of UIDs
nuccore	binary	asn.1	none	binary ASN.1
nuccore	fasta	text	fasta	FASTA
nuccore	fasta	text	fasta_cds_aa	CDS protein FASTA
nuccore	fasta	text	fasta_cds_na	CDS nucleotide FASTA
nuccore	genbank	text	gb	GenBank flat file
nuccore	genbank	text	gbwithparts	GenBank flat file with full sequence (contigs)
nuccore	tabular	text	acc	Accession number(s)
nuccore	txt	text	ft	Feature table
nuccore	tabular	text	seqid	SeqID string
nuccore	tabular	text	uilist	List of UIDs
nuccore	txt	text	none	text ASN.1
nuccore	xml	xml	docsum	Document summary
nuccore	xml	xml	fasta	TinySeq
nuccore	xml	xml	gb	GBSeq
nuccore	xml	xml	gbc	INSDSeq
nuccore	xml	xml	native	Full record
nuccore	xml	xml	uilist	List of UIDs
nucest	binary	asn.1	none	binary ASN.1
nucest	fasta	text	fasta	FASTA
nucest	genbank	text	gb	GenBank flat file
nucest	tabular	text	acc	Accession number(s)
nucest	tabular	text	seqid	SeqID string
nucest	tabular	text	uilist	List of UIDs
nucest	txt	text	est	EST report
nucest	txt	text	none	text ASN.1
nucest	xml	xml	docsum	Document summary
nucest	xml	xml	fasta	TinySeq
nucest	xml	xml	gb	GBSeq
nucest	xml	xml	gbc	INSDSeq
nucest	xml	xml	native	Full record
nucest	xml	xml	uilist	List of UIDs
nucgss	binary	asn.1	none	binary ASN.1
nucgss	fasta	text	fasta	FASTA
nucgss	genbank	text	gb	GenBank flat file
nucgss	tabular	text	acc	Accession number(s)
nucgss	tabular	text	seqid	SeqID string
nucgss	tabular	text	uilist	List of UIDs
nucgss	txt	text	gss	GSS report
nucgss	txt	text	none	text ASN.1
nucgss	xml	xml	docsum	Document summary
nucgss	xml	xml	fasta	TinySeq
nucgss	xml	xml	gb	GBSeq
nucgss	xml	xml	gbc	INSDSeq
nucgss	xml	xml	native	Full record
nucgss	xml	xml	uilist	List of UIDs
pmc	tabular	text	uilist	List of UIDs
pmc	txt	text	medline	MEDLINE
pmc	xml	xml	docsum	Document summary
pmc	xml	xml	none	FULL
pmc	xml	xml	uilist	List of UIDs
popset	binary	asn.1	none	binary ASN.1
popset	fasta	text	fasta	FASTA
popset	genbank	text	gb	GenBank flat file
popset	tabular	text	acc	Accession number(s)
popset	tabular	text	seqid	SeqID string
popset	tabular	text	uilist	List of UIDs
popset	txt	text	none	text ASN.1
popset	xml	xml	docsum	Document summary
popset	xml	xml	fasta	TinySeq
popset	xml	xml	gb	GBSeq
popset	xml	xml	gbc	INSDSeq
popset	xml	xml	native	Full record
popset	xml	xml	uilist	List of UIDs
protein	binary	asn.1	none	binary ASN.1
protein	fasta	text	fasta	FASTA
protein	tabular	text	acc	Accession number(s)
protein	txt	text	ft	Feature table
protein	tabular	text	seqid	SeqID string
protein	tabular	text	uilist	List of UIDs
protein	txt	text	gp	GenPept flat file
protein	txt	text	none	text ASN.1
protein	xml	xml	docsum	Document summary
protein	xml	xml	fasta	TinySeq
protein	xml	xml	gp	GBSeq
protein	xml	xml	gpc	INSDSeq
protein	xml	xml	ipg	Identical Protein
protein	xml	xml	native	Full record
protein	xml	xml	uilist	List of UIDs
pubmed	tabular	text	uilist	List of UIDs
pubmed	txt	asn.1	none	text ASN.1
pubmed	txt	text	abstract	Abstract
pubmed	txt	text	medline	MEDLINE
pubmed	xml	xml	docsum	Document summary
pubmed	xml	xml	none	Full
pubmed	xml	xml	uilist	List of UIDs
sequences	fasta	text	fasta	FASTA
sequences	tabular	text	acc	Accession number(s)
sequences	tabular	text	seqid	SeqID string
sequences	tabular	text	uilist	List of UIDs
sequences	txt	text	none	text ASN.1
sequences	xml	xml	docsum	Document summary
sequences	xml	xml	uilist	List of UIDs
sequences	xml	none	none	Full
snp	fasta	text	fasta	FASTA
snp	json	json	docsum	Document summary
snp	json	json	uilist	List of UIDs
snp	tabular	text	ssexemplar	SS Exemplar list
snp	tabular	text	uilist	List of UIDs
snp	txt	asn.1	none	text ASN.1
snp	txt	text	chr	Chromosome report
snp	txt	text	docset	Summary
snp	txt	text	flt	Flat file
snp	txt	text	rsr	RS Cluster report
snp	xml	xml	docsum	Document summary
snp	xml	xml	none	XML
snp	xml	xml	uilist	List of UIDs
sra	tabular	text	uilist	List of UIDs
sra	xml	xml	docsum	Document summary
sra	xml	xml	full	Full
taxonomy	tabular	text	uilist	List of UIDs
taxonomy	xml	xml	none	Full
taxonomy	xml	xml	docsum	Document summary
taxonomy	xml	xml	uilist	List of UIDs
