== RSeQC Galaxy Wrapper ==

This is a Galaxy wrapper for the RSeQC RNA-Seq QC package.

** Installation **

Installation from a tool shed provides the necessary tool dependencies, R, numpy, and RSeQC.

Otherwise, make sure that R and the RSeQC scripts are in the path and run under the Galaxy environment.
Move the xml files to a subdirectory of your tools directory and add lines in tool_conf.xml to point to them.
Restart the Galaxy server.

Requires Python 2.7

** Attribution **

The RSeQC package and associated documentation can be found at: http://rseqc.sourceforge.net/

The galaxy wrapper code was written by
    Nilesh Kavthekar, School of Engineering and Applied Sciences, University of Pennsylvania, Class of 2016
Modified by
    Lance Parsons, Lewis-Sigler Institute for Integrative Genomics, Princeton University,
    Bjorn Gruning, University of Freiburg, bjoern.gruening@gmail.com
