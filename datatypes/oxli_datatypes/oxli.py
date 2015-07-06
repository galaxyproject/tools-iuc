"""
k-mer count and presence
"""

from galaxy.datatypes.binary import Binary
import binascii

import logging

log = logging.getLogger(__name__)


class OxliBinary(Binary):

    def __init__(self, **kwd):
        OxliBinary.__init__(self, **kwd)

    def sniff(self, filename, filetype):
        try:
            with open(filename) as fileobj:
                header = fileobj.read(4)
                if binascii.b2a_hex(header) == binascii.hexlify('OXLI'):
                    fileobj.seek(1)
                    ftype = fileobj.read(1)
                    if binascii.b2a_hex(ftype) == filetype:
                        return True
            return False
        except IOError:
            return False


class Count(OxliBinary):

    def __init__(self, **kwd):
        OxliBinary.__init__(self, **kwd)

    def sniff(self, filename):
        return OxliBinary.sniff(self, filename, "01")


class Presence(OxliBinary):

    def __init__(self, **kwd):
        OxliBinary.__init__(self, **kwd)

    def sniff(self, filename):
        return OxliBinary.sniff(self, filename, "02")


Binary.register_sniffable_binary_ext("ct", "ct", Count)
Binary.register_sniffable_binary_ext("pt", "pt", Presence)
