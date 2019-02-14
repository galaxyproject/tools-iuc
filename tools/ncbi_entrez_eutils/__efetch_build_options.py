#!/usr/bin/env python
# Daniel Blankenberg
# Creates the options for tool interface
from __future__ import print_function

import re

# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi
db_list = '''
<DbName>annotinfo</DbName>
<DbName>assembly</DbName>
<DbName>bioproject</DbName>
<DbName>biosample</DbName>
<DbName>biosystems</DbName>
<DbName>blastdbinfo</DbName>
<DbName>books</DbName>
<DbName>cdd</DbName>
<DbName>clinvar</DbName>
<DbName>clone</DbName>
<DbName>dbvar</DbName>
<DbName>gap</DbName>
<DbName>gapplus</DbName>
<DbName>gds</DbName>
<DbName>gencoll</DbName>
<DbName>gene</DbName>
<DbName>genome</DbName>
<DbName>geoprofiles</DbName>
<DbName>grasp</DbName>
<DbName>gtr</DbName>
<DbName>homologene</DbName>
<DbName>medgen</DbName>
<DbName>mesh</DbName>
<DbName>ncbisearch</DbName>
<DbName>nlmcatalog</DbName>
<DbName>nuccore</DbName>
<DbName>nucest</DbName>
<DbName>nucgss</DbName>
<DbName>nucleotide</DbName>
<DbName>omim</DbName>
<DbName>orgtrack</DbName>
<DbName>pcassay</DbName>
<DbName>pccompound</DbName>
<DbName>pcsubstance</DbName>
<DbName>pmc</DbName>
<DbName>popset</DbName>
<DbName>probe</DbName>
<DbName>protein</DbName>
<DbName>proteinclusters</DbName>
<DbName>pubmed</DbName>
<DbName>pubmedhealth</DbName>
<DbName>seqannot</DbName>
<DbName>snp</DbName>
<DbName>sra</DbName>
<DbName>structure</DbName>
<DbName>taxonomy</DbName>
<DbName>unigene</DbName>'''.replace( "<DbName>", "").replace( "</DbName>", "").split("\n")


help = '''  (all)
                 docsum             xml      Document Summary
                 docsum             json     Document Summary
                 full               text     Full Document
                 uilist             xml      Unique Identifier List
                 uilist             text     Unique Identifier List
                 full               xml      Full Document

  bioproject
                 native                      BioProject Report
                 native             xml      RecordSet

  biosample
                 native                      BioSample Report
                 native             xml      BioSampleSet

  biosystems
                 native             xml      Sys-set

  gds
                 native             xml      RecordSet
                 summary            text     Summary

  gene
                 gene_table         xml      Gene Table
                 native             text     Gene Report
                 native             asn.1    Entrezgene
                 native             xml      Entrezgene-Set
                 tabular            tabular  Tabular Report

  homologene
                 alignmentscores    text     Alignment Scores
                 fasta              fasta    FASTA
                 homologene         text     Homologene Report
                 native             text     Homologene List
                 native             asn.1    HG-Entry
                 native             xml      Entrez-Homologene-Set

  mesh
                 full               text     Full Record
                 native             text     MeSH Report
                 native             xml      RecordSet

  nlmcatalog
                 native             text     Full Record
                 native             xml      NLMCatalogRecordSet

  pmc
                 medline            text     MEDLINE
                 native             xml      pmc-articleset

  pubmed
                 abstract           xml      Abstract
                 medline            text     MEDLINE
                 native             asn.1    Pubmed-entry
                 native             xml      PubmedArticleSet

  (sequences)
                 acc                text     Accession Number
                 est                xml      EST Report
                 fasta              fasta    FASTA
                 fasta              xml      TinySeq
                 fasta_cds_aa       fasta    CDS Products
                 fasta_cds_na       fasta    Coding Regions
                 ft                 text     Feature Table
                 gb                 text     GenBank Flatfile
                 gb                 xml      GBSet
                 gbc                xml      INSDSet
                 gbwithparts        text     GenBank with Contig Sequences
                 gene_fasta         fasta    FASTA of Gene
                 gp                 text     GenPept Flatfile
                 gp                 xml      GBSet
                 gpc                xml      INSDSet
                 gss                text     GSS Report
                 ipg                text     Identical Protein Report
                 ipg                xml      IPGReportSet
                 native             text     Seq-entry
                 native             xml      Bioseq-set
                 seqid              asn.1    Seq-id

  snp
                 chr                text     Chromosome Report
                 docset             text     Summary
                 fasta              fasta    FASTA
                 flt                text     Flat File
                 native             asn.1    Rs
                 native             xml      ExchangeSet
                 rsr                tabular  RS Cluster Report
                 ssexemplar         text     SS Exemplar List

  sra
                 native             xml      EXPERIMENT_PACKAGE_SET
                 runinfo            xml      SraRunInfo

  structure
                 mmdb               asn.1    Ncbi-mime-asn1 strucseq
                 native             text     MMDB Report
                 native             xml      RecordSet

  taxonomy
                 native             text     Taxonomy List
                 native             xml      TaxaSet'''.split("\n")


db = {}
for db_name in db_list:
    db[db_name] = []

section = None
for line in help:
    line = re.split(r'\s{2,}', line.strip())
    # Ignore empties
    if len(line) == 0:
        continue
    # Section headers have one item
    elif len(line) == 1:
        section = line[0]
        db[section] = []
    # Format lines have 2+
    elif len(line) == 2:
        parent_format = line[0]
        description = line[1]

        if parent_format not in db[section]:
            db[section].append((parent_format, None, description))
    elif len(line) == 3:
        parent_format = line[0]
        format_modifier = line[1]
        description = line[2]

        if parent_format not in db[section]:
            db[section].append((parent_format, format_modifier, description))


all_formats = db['(all)']
del db['(all)']
sequences_formats = db['(sequences)']
del db['(sequences)']
del db['']

for key in db:
    db[key] += all_formats

for key in ('nuccore', 'nucest', 'nucgss', 'nucleotide'):
    db[key] += sequences_formats

MACRO_TPL = '''

'''

WHEN_TPL = '''      <when value="{format}">
        <param name="output_format" type="select" label="Output Format">
          {format_options}
        </param>
      </when>'''

FORMAT_OPTION_TPL = '''<option value="{name_type}">{name_type_human}</option>'''

format_names = {}

print('''  <xml name="db">
    <conditional name="db">
      <expand macro="dbselect" />''')
for key in sorted(db):
    format_options = []

    for (parent_format, format_modifier, description) in sorted(db[key]):
        name_human = description
        if format_modifier:
            name_human += ' (%s)' % format_modifier
        format_string = '%s-%s' % (parent_format, format_modifier)

        format_options.append(FORMAT_OPTION_TPL.format(
            name_type=format_string,
            name_type_human=name_human,
        ))

        format_names[format_string] = format_modifier

    print(WHEN_TPL.format(
        format=key,
        format_options='\n          '.join(format_options)
    ))

print('''    </conditional>
  </xml>''')

CHANGE_FORMAT_TPL = '''
  <xml name="efetch_formats">
    <change_format>
      {formats}
    </change_format>
  </xml>
'''

CHANGE_FORMAT_WHEN_TPL = '''<when input="output_format" value="{key}" format="{value}"/>'''
# Format options


whens = []
for (k, v) in format_names.items():
    if v is None:
        v = 'text'
    elif v == 'asn.1':
        v = 'asn1'

    whens.append(CHANGE_FORMAT_WHEN_TPL.format(
        key=k, value=v
    ))

print(CHANGE_FORMAT_TPL.format(formats='\n      '.join(whens)))
