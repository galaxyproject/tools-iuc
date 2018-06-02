"""
SnpSift dbNSFP datatypes
"""
import gzip
import logging
import os
import os.path
import sys
import traceback

from galaxy.datatypes.data import Text
from galaxy.datatypes.metadata import MetadataElement

log = logging.getLogger(__name__)


class SnpSiftDbNSFP( Text ):
    """Class describing a dbNSFP database prepared fpr use by SnpSift dbnsfp """
    MetadataElement( name='reference_name', default='dbSNFP', desc='Reference Name', readonly=True, visible=True, set_in_upload=True, no_value='dbSNFP' )
    MetadataElement( name="bgzip", default=None, desc="dbNSFP bgzip", readonly=True, visible=True, no_value=None )
    MetadataElement( name="index", default=None, desc="Tabix Index File", readonly=True, visible=True, no_value=None)
    MetadataElement( name="annotation", default=[], desc="Annotation Names", readonly=True, visible=True, no_value=[] )
    file_ext = "snpsiftdbnsfp"
    composite_type = 'auto_primary_file'
    allow_datatype_change = False
    """
    ## The dbNSFP file is a tabular file with 1 header line
    ## The first 4 columns are required to be: chrom	pos	ref	alt
    ## These match columns 1,2,4,5 of the VCF file
    ## SnpSift requires the file to be block-gzipped and the indexed with samtools tabix
    ## Example:
    ## Compress using block-gzip algorithm
    bgzip dbNSFP2.3.txt
    ## Create tabix index
    tabix -s 1 -b 2 -e 2 dbNSFP2.3.txt.gz
    """
    def __init__( self, **kwd ):
        Text.__init__( self, **kwd )
        self.add_composite_file('%s.grp', description='Group File', substitute_name_with_metadata='reference_name', is_binary=False)
        self.add_composite_file('%s.ti', description='', substitute_name_with_metadata='reference_name', is_binary=False)

    def init_meta( self, dataset, copy_from=None ):
        Text.init_meta( self, dataset, copy_from=copy_from )

    def generate_primary_file(self, dataset=None):
        """
        This is called only at upload to write the html file
        cannot rename the datasets here - they come with the default unfortunately
        """
        self.regenerate_primary_file(dataset)

    def regenerate_primary_file(self, dataset):
        """
        cannot do this until we are setting metadata
        """
        annotations = "dbNSFP Annotations: %s\n" % ','.join(dataset.metadata.annotation)
        f = open(dataset.file_name, 'a')
        if dataset.metadata.bgzip:
            bn = dataset.metadata.bgzip
            f.write(bn)
            f.write('\n')
        f.write(annotations)
        f.close()

    def set_meta( self, dataset, overwrite=True, **kwd ):
        try:
            efp = dataset.extra_files_path
            if os.path.exists(efp):
                flist = os.listdir(efp)
                for i, fname in enumerate(flist):
                    if fname.endswith('.gz'):
                        dataset.metadata.bgzip = fname
                        try:
                            fh = gzip.open(os.path.join(efp, fname), 'r')
                            buf = fh.read(5000)
                            lines = buf.splitlines()
                            headers = lines[0].split('\t')
                            dataset.metadata.annotation = headers[4:]
                        except Exception as e:
                            log.warn("set_meta fname: %s  %s" % (fname, str(e)))
                            traceback.print_stack(file=sys.stderr)
                        finally:
                            fh.close()
                    if fname.endswith('.tbi'):
                        dataset.metadata.index = fname
            self.regenerate_primary_file(dataset)
        except Exception as e:
            log.warn("set_meta fname: %s  %s" % (dataset.file_name if dataset and dataset.file_name else 'Unkwown', str(e)))
            traceback.print_stack(file=sys.stderr)


if __name__ == '__main__':
    import doctest
    doctest.testmod(sys.modules[__name__])
