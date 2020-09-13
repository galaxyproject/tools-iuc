The Galaxy tool wrappers contained in this tool shed repository rely on software developed by
the NCBI: https://github.com/ncbi/sra-tools.

NCBI Sequence Read Archive Toolkit: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software

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
