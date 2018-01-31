"""
SnpEff datatypes
"""
import gzip
import logging
import os
import os.path
import re

from galaxy.datatypes.data import Text
from galaxy.datatypes.metadata import MetadataElement

log = logging.getLogger(__name__)


class SnpEffDb( Text ):
    """Class describing a SnpEff genome build"""
    file_ext = "snpeffdb"
    MetadataElement( name="genome_version", default=None, desc="Genome Version", readonly=True, visible=True, no_value=None )
    MetadataElement( name="snpeff_version", default="SnpEff4.0", desc="SnpEff Version", readonly=True, visible=True, no_value=None )
    MetadataElement( name="regulation", default=[], desc="Regulation Names", readonly=True, visible=True, no_value=[], optional=True)
    MetadataElement( name="annotation", default=[], desc="Annotation Names", readonly=True, visible=True, no_value=[], optional=True)

    def __init__( self, **kwd ):
        Text.__init__( self, **kwd )

    # The SnpEff version line was added in SnpEff version 4.1
    def getSnpeffVersionFromFile(self, path):
        snpeff_version = None
        try:
            fh = gzip.open(path, 'rb')
            buf = fh.read(100)
            lines = buf.splitlines()
            m = re.match('^(SnpEff)\s+(\d+\.\d+).*$', lines[0].strip())
            if m:
                snpeff_version = m.groups()[0] + m.groups()[1]
            fh.close()
        except Exception:
            pass
        return snpeff_version

    def set_meta( self, dataset, **kwd ):
        Text.set_meta(self, dataset, **kwd )
        data_dir = dataset.extra_files_path
        # search data_dir/genome_version for files
        regulation_pattern = 'regulation_(.+).bin'
        #  annotation files that are included in snpEff by a flag
        annotations_dict = {'nextProt.bin': '-nextprot', 'motif.bin': '-motif'}
        regulations = []
        annotations = []
        genome_version = None
        snpeff_version = None
        if data_dir and os.path.isdir(data_dir):
            for root, dirs, files in os.walk(data_dir):
                for fname in files:
                    if fname.startswith('snpEffectPredictor'):
                        # if snpEffectPredictor.bin download succeeded
                        genome_version = os.path.basename(root)
                        dataset.metadata.genome_version = genome_version
                    else:
                        m = re.match(regulation_pattern, fname)
                        if m:
                            name = m.groups()[0]
                            regulations.append(name)
                        elif fname in annotations_dict:
                            value = annotations_dict[fname]
                            name = value.lstrip('-')
                            annotations.append(name)
            dataset.metadata.regulation = regulations
            dataset.metadata.annotation = annotations
            try:
                fh = file(dataset.file_name, 'w')
                fh.write("%s\n" % genome_version if genome_version else 'Genome unknown')
                fh.write("%s\n" % snpeff_version if snpeff_version else 'SnpEff version unknown')
                if annotations:
                    fh.write("annotations: %s\n" % ','.join(annotations))
                if regulations:
                    fh.write("regulations: %s\n" % ','.join(regulations))
                fh.close()
            except Exception:
                pass
