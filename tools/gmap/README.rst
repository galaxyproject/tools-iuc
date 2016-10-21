Galaxy Wrapper for GMAP
=======================

This wrapper is copyright 2011, 2013 by Jim Johnson (Minnesota Supercomputing
Institute, University of Minnesota).
Revisions copyright 2016 by Peter Cock (The James Hutton Institute, UK).

This Galaxy wrapper is available from the Galaxy Tool Shed at:

toolshed.g2.bx.psu.edu/view/jjohnson/gmap

GMAP applications and citation info are available from:
ttp://research-pub.gene.com/gmap/

Manual GMAP installation instructions are in the README file in the download,
and online: http://research-pub.gene.com/gmap/src/README


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v2.0.0  - Initial release, for GMAP version 2011-10-07 (by JJ, October 2011).
v2.0.1  - Updated for GMAP version 2011-11-30 (by JJ, November 2011).
v3.0.0  - Work in progress for GMAP version 2013-05-09 (by JJ, June 2013),
          never released to the main Galaxy Tool Shed.
v3.0.1  - Various style updates (Peter Cock, October 2016).
======= ======================================================================


Automated Installation
======================

This should be straightforward using the Galaxy Tool Shed, which should be
able to automatically install GMAP as well.



Manual Installation
===================

GMAP and  GSNAP use added datatypes:

 -  add datatype definition file: ``lib/galaxy/datatypes/gmap.py``

 -  add the following import line to:  ``lib/galaxy/datatypes/registry.py``::
  
        import gmap # added for gmap tools

 -  add to ``datatypes_conf.xml``::

        <!-- Start GMAP Datatypes -->
        <datatype extension="gmapdb" type="galaxy.datatypes.gmap:GmapDB"  display_in_upload="False"/>
        <datatype extension="gmapsnpindex" type="galaxy.datatypes.gmap:GmapSnpIndex"  display_in_upload="False"/>
        <datatype extension="iit" type="galaxy.datatypes.gmap:IntervalIndexTree"  display_in_upload="True"/>
        <datatype extension="splicesites.iit" type="galaxy.datatypes.gmap:SpliceSitesIntervalIndexTree"  display_in_upload="True"/>
        <datatype extension="introns.iit" type="galaxy.datatypes.gmap:IntronsIntervalIndexTree"  display_in_upload="True"/>
        <datatype extension="snps.iit" type="galaxy.datatypes.gmap:SNPsIntervalIndexTree"  display_in_upload="True"/>
        <datatype extension="tally.iit" type="galaxy.datatypes.gmap:TallyIntervalIndexTree"  display_in_upload="True"/>
        <datatype extension="gmap_annotation" type="galaxy.datatypes.gmap:IntervalAnnotation"  display_in_upload="False"/>
        <datatype extension="gmap_splicesites" type="galaxy.datatypes.gmap:SpliceSiteAnnotation"  display_in_upload="True"/>
        <datatype extension="gmap_introns" type="galaxy.datatypes.gmap:IntronAnnotation"  display_in_upload="True"/>
        <datatype extension="gmap_snps" type="galaxy.datatypes.gmap:SNPAnnotation"  display_in_upload="True"/>
        <datatype extension="gsnap_tally" type="galaxy.datatypes.gmap:TallyAnnotation"  display_in_upload="True"/>
        <datatype extension="gsnap" type="galaxy.datatypes.gmap:GsnapResult"  display_in_upload="True"/>
        <!-- End GMAP Datatypes -->

Tools:
 - GMAP_Build - create a GmapDB set of index files for a reference sequence and optional set of annotations
 - GMAP - map sequences to a reference sequence GmapDB index
 - GSNAP - align sequences to a reference and detect splicing

Add to ``tool_conf.xml`` (probably in the "NGS: Mapping" section)::

   <tool file="gmap/gmap.xml" />
   <tool file="gmap/gsnap.xml" />
   <tool file="gmap/gmap_build.xml" />
   <tool file="gmap/snpindex.xml" />
   <tool file="gmap/iit_store.xml" />

Admin built cached gmapdb indexes defined in ``tool-data/gmap_indices.loc``


TODO
====

 - Add classes to ``gmap.py``

   - CmetIndex - an index created by cmetindex
   - AtoiIndex - an index created by atoiindex

 - Add tally creation

   - gsnap default output -> gsnap_tally -> iit_store

 - Add goby support

   - Should add separate tools and datatypes for goby
   - GSNAP goby output relies on goby input, might be better to have a separate gsnap tool for goby

 - Possibly add Tools:

   - get_genome - retrieves from a gmapdb
   - cmetindex - create methylcytosine index
   - atoiindex - create  A-to-I RNA editing index
