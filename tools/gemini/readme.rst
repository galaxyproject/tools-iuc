=========================
Galaxy wrapper for GEMINI
=========================


GEMINI: a flexible framework for exploring genome variation

GEMINI (GEnome MINIng) is designed to be a flexible framework for exploring genetic variation in the context of 
the wealth of genome annotations available for the human genome. By placing genetic variants, sample genotypes, 
and useful genome annotations into an integrated database framework, GEMINI provides a simple, flexible, yet very 
powerful system for exploring genetic variation for for disease and population genetics.

Using the GEMINI framework begins by loading a VCF file into a database. Each variant is automatically 
annotated by comparing it to several genome annotations from source such as ENCODE tracks, UCSC tracks, 
OMIM, dbSNP, KEGG, and HPRD. All of this information is stored in portable SQLite database that allows 
one to explore and interpret both coding and non-coding variation using “off-the-shelf” tools or an 
enhanced SQL engine.

Please also see the original [manuscript](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003153).


============
Installation
============

It is recommended to install this wrapper via the `Galaxy Tool Shed`.

.. _`Galaxy Tool Shed`:  https://testtoolshed.g2.bx.psu.edu/view/iuc/gemini


=======
History
=======
- 0.9.1: Initial public release


====================
Detailed description
====================

View the original GEMINI documentation: http://gemini.readthedocs.org/en/latest/index.html


===============================
Wrapper Licence (MIT/BSD style)
===============================

Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the
names of the contributors or copyright holders not be used in
advertising or publicity pertaining to distribution of the software
without specific prior permission.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.

