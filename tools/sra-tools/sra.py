"""
NCBI sra class
"""
import logging
import binascii
from galaxy.datatypes.data import nice_size
from galaxy.datatypes.binary import Binary

log = logging.getLogger(__name__)

class Sra(Binary):
    """ Sequence Read Archive (SRA) """
    file_ext = 'sra'

    def __init__( self, **kwd ):
        Binary.__init__( self, **kwd )

    def sniff( self, filename ):
        """ The first 8 bytes of any NCBI sra file is 'NCBI.sra', and the file is binary.
        For details about the format, see http://www.ncbi.nlm.nih.gov/books/n/helpsra/SRA_Overview_BK/#SRA_Overview_BK.4_SRA_Data_Structure
        """
        try:
            header = open(filename).read(8)
            if binascii.b2a_hex(header) == binascii.hexlify('NCBI.sra'):
                return True
            else:
                return False
        except:
            return False

    def set_peek(self, dataset, is_multi_byte=False):
        if not dataset.dataset.purged:
            dataset.peek  = 'Binary sra file'
            dataset.blurb = nice_size(dataset.get_size())
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'

    def display_peek(self, dataset):
        try:
            return dataset.peek
        except:
            return 'Binary sra file (%s)' % (nice_size(dataset.get_size()))
            
Binary.register_sniffable_binary_format('sra', 'sra', Sra)
