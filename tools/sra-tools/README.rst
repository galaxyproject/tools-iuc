The Galaxy tool wrappers contained in this tool shed repository rely on software developed by
the NCBI: http://github.com/ncbi/sra-tools.

NCBI Sequence Read Archive: http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software.
Use of SRA Toolkit software herein should comply with the GPL v2 or greater.

Copyright (C) 2013  Matthew Shirley

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# INSTALLATION

This software release was designed to install using the Galaxy toolshed under Linux and MacOS operating systems on Intel x86-compatible 32/64 bit architectures.

*Build Requirements*

- bash
- make
- gcc
- g++
- libxml2

On a Debian OS use:

    apt-get install build-essential libxml2-dev

On a Mac with [command line tools](https://developer.apple.com/downloads/index.action) installed:

    brew install libxml2

# Installation of Aspera connect ascp binary

The sra-tools suite is ready to benefit from increased transfer speed and reliability by using Aspera Connect ascp.
To benefit, download the ascp commandline client, and place ascp and the required ssh keys into a PATH accessible to galaxys job handler.

A convenience package for linux and OS X is available at https://toolshed.g2.bx.psu.edu/view/mvdbeek/package_ascp_3/e109f0ec22c3 .
It suffices to copy the contents of the $INSTALL_DIR/bin into PATH.

Alternatively go to http://downloads.asperasoft.com/connect2/ .

Aspera connect is not provided by the IUC due to its closed-source nature.

# Firewall settings for highspeed transfer

To benefit from increased transfer speeds using ascp3 your local firewall must permit UDP data transfer in both
directions on ports 33001-33009 for the following IP ranges:

    130.14.*.*

    165.112.*.*

The firewall must also allow ssh traffic outbound to NCBI.
The wrapper will fall back to http download if these requirements are not met.

CONTROLLED-ACCESS DATA

Encrypted, controlled-access data is not supported.
