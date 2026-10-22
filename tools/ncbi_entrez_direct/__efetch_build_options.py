#!/usr/bin/env python

# Daniel Blankenberg
# Creates the options for tool interface

# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi
db_list = '''<DbName>pubmed</DbName>
<DbName>protein</DbName>
<DbName>nuccore</DbName>
<DbName>nucleotide</DbName>
<DbName>nucgss</DbName>
<DbName>nucest</DbName>
<DbName>structure</DbName>
<DbName>genome</DbName>
<DbName>annotinfo</DbName>
<DbName>assembly</DbName>
<DbName>bioproject</DbName>
<DbName>biosample</DbName>
<DbName>blastdbinfo</DbName>
<DbName>books</DbName>
<DbName>cdd</DbName>
<DbName>clinvar</DbName>
<DbName>clone</DbName>
<DbName>gap</DbName>
<DbName>gapplus</DbName>
<DbName>grasp</DbName>
<DbName>dbvar</DbName>
<DbName>gene</DbName>
<DbName>gds</DbName>
<DbName>geoprofiles</DbName>
<DbName>homologene</DbName>
<DbName>medgen</DbName>
<DbName>mesh</DbName>
<DbName>ncbisearch</DbName>
<DbName>nlmcatalog</DbName>
<DbName>omim</DbName>
<DbName>orgtrack</DbName>
<DbName>pmc</DbName>
<DbName>popset</DbName>
<DbName>probe</DbName>
<DbName>proteinclusters</DbName>
<DbName>pcassay</DbName>
<DbName>biosystems</DbName>
<DbName>pccompound</DbName>
<DbName>pcsubstance</DbName>
<DbName>pubmedhealth</DbName>
<DbName>seqannot</DbName>
<DbName>snp</DbName>
<DbName>sra</DbName>
<DbName>taxonomy</DbName>
<DbName>unigene</DbName>
<DbName>gencoll</DbName>
<DbName>gtr</DbName>'''.replace("<DbName>", "").replace("</DbName>", "").split("\n")


help = '''  (all)
                 docsum                      DocumentSummarySet XML
                 docsum             json     DocumentSummarySet JSON
                 full                        Same as native except for mesh
                 uid                         Unique Identifier List
                 url                         Entrez URL
                 xml                         Same as -format full -mode xml

  bioproject
                 native                      BioProject Report
                 native             xml      RecordSet XML

  biosample
                 native                      BioSample Report
                 native             xml      BioSampleSet XML

  biosystems
                 native             xml      Sys-set XML

  gds
                 native             xml      RecordSet XML
                 summary                     Summary

  gene
                 gene_table                  Gene Table
                 native                      Gene Report
                 native             asn.1    Entrezgene ASN.1
                 native             xml      Entrezgene-Set XML
                 tabular                     Tabular Report

  homologene
                 alignmentscores             Alignment Scores
                 fasta                       FASTA
                 homologene                  Homologene Report
                 native                      Homologene List
                 native             asn.1    HG-Entry ASN.1
                 native             xml      Entrez-Homologene-Set XML

  mesh
                 full                        Full Record
                 native                      MeSH Report
                 native             xml      RecordSet XML

  nlmcatalog
                 native                      Full Record
                 native             xml      NLMCatalogRecordSet XML

  pmc
                 medline                     MEDLINE
                 native             xml      pmc-articleset XML

  pubmed
                 abstract                    Abstract
                 medline                     MEDLINE
                 native             asn.1    Pubmed-entry ASN.1
                 native             xml      PubmedArticleSet XML

  (sequences)
                 acc                         Accession Number
                 est                         EST Report
                 fasta                       FASTA
                 fasta              xml      TinySeq XML
                 fasta_cds_aa                FASTA of CDS Products
                 fasta_cds_na                FASTA of Coding Regions
                 ft                          Feature Table
                 gb                          GenBank Flatfile
                 gb                 xml      GBSet XML
                 gbc                xml      INSDSet XML
                 gbwithparts                 GenBank with Contig Sequences
                 gene_fasta                  FASTA of Gene
                 gp                          GenPept Flatfile
                 gp                 xml      GBSet XML
                 gpc                xml      INSDSet XML
                 gss                         GSS Report
                 ipg                         Identical Protein Report
                 ipg                xml      IPGReportSet XML
                 native             text     Seq-entry ASN.1
                 native             xml      Bioseq-set XML
                 seqid                       Seq-id ASN.1

  snp
                 chr                         Chromosome Report
                 docset                      Summary
                 fasta                       FASTA
                 flt                         Flat File
                 native             asn.1    Rs ASN.1
                 native             xml      ExchangeSet XML
                 rsr                         RS Cluster Report
                 ssexemplar                  SS Exemplar List

  sra
                 native             xml      EXPERIMENT_PACKAGE_SET XML
                 runinfo            xml      SraRunInfo XML

  structure
                 mmdb                        Ncbi-mime-asn1 strucseq ASN.1
                 native                      MMDB Report
                 native             xml      RecordSet XML

  taxonomy
                 native                      Taxonomy List
                 native             xml      TaxaSet XML'''.split("\n")

db = {}
name = None
all = "(all)"
for line in help:
    if line.strip() and line[2] != ' ':
        name = line.strip()
        db[name] = {}
    elif line.strip():
        format = line[0:len("                 docsum             ")].strip()
        mode = line[len("                 docsum             "):len("                 docsum             json     ")].strip()
        if format not in db[name]:
            db[name][format] = []
        db[name][format].append(mode)

for name in db_list:
    if name not in db:
        db[name] = {}

db["sequences"] = db["(sequences)"]
del db["(sequences)"]

print('<conditional name="db">')
print('    <param name="db" type="select" label="Database" argument="-db">')
for name in sorted(db.keys()):
    if name == all:
        continue
    print('        <option value="%s">%s</option>' % (name, name))
print('        <option value="">Manual Entry</option>')
print('    </param>')

for name in sorted(db.keys()):
    if name == all:
        continue
    my_dict = db[all].copy()

    for format, modes in db[name].items():
        if format in my_dict:
            for mode in modes:
                if mode not in my_dict[format]:
                    my_dict[format].append(mode)
        else:
            my_dict[format] = modes
    if "" not in my_dict:
        my_dict[""] = [""]
    print('    <when value="%s">' % name)
    print('        <conditional name="format">')
    print('            <param name="format" type="select" label="Format" argument="-format">')
    for format in sorted(my_dict.keys()):
        print('                <option value="%s">%s</option>' % (format, format or "None"))
    print('            </param>')
    for format in sorted(my_dict.keys()):
        print('            <when value="%s">' % format)
        print('                <param name="mode" type="select" label="Mode" argument="-mode">')
        if "" not in my_dict[format]:
            my_dict[format].append("")
        for mode in sorted(my_dict[format]):
            print('                    <option value="%s">%s</option>' % (mode, mode or "None"))
        print('                </param>')
        print('            </when>')
    print('        </conditional>')
    print('    </when>')
print('    <when value="">')
print('        <param name="db_manual" type="text" label="Database" argument="-db"/>')
print('        <param name="format" type="text" label="Format" argument="-format"/>')
print('        <param name="mode" type="text" label="Mode" argument="-mode"/>')
print('    </when>')
print('</conditional>')
