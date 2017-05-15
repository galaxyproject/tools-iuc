#!/usr/bin/env python
'''
Extract the blocks from a BGZip file.

BAM files are stored as blocks in a bgzip archive. This class
will load the bgzip archive and output the block information.
'''

import os
import struct
import sys


class BGZip(object):
    def __init__(self, fname):
        self.fname = fname
        self.pos = 0
        self.fileobj = open(self.fname)
        self.fsize = os.stat(self.fname).st_size
        self.cur_chunk = 0

        self.cpos = 0
        self.cdata = 0

    def close(self):
        self.fileobj.close()

    def seek(self, offset):
        bgz_offset = offset >> 16
        chunk_offset = offset & 0xFFFF

        self.fileobj.seek(bgz_offset)
        self.read_chunk()
        self.chunk_pos = chunk_offset

    def read(self, amount, whence=0):
        if whence not in [0, 1]:
            print "Bad Whence value!: %s" % whence
            return

        if whence == 0:
            self.seek(0, 0)

        # read into chunk, if not enough data in chunk, read next chunk
        ret = ''
        while amount and self.pos < self.fsize:
            if len(self.cdata) - self.cpos < amount:
                ret += self.cdata[self.cpos:self.cpos + amount]
                self.cpos += amount
                return ret
            else:
                ret += self.cdata[self.cpos:]
                amount = amount - len(ret)
                self.read_chunk()
        return ret

    def read_chunk(self):
        self.fileobj.seek(10)
        id1, id2, cm, flg, mtime, xfl, os, xlen = self._read_fields('<BBBBIBBH')
        subpos = 0
        bsize = 0

        while subpos < xlen:
            si1, si2, slen = self._read_fields('<BBH')
            if si1 == 66 and si2 == 67:
                bsize, = self._read_fields('<H')
            else:
                self.fileobj.seek(slen, 1)
                self.pos += slen

            subpos += 6 + slen

        cdata_size = bsize - xlen - 19
        self.cdata = self.fileobj.read(cdata_size)  # inflate value
        self.fileobj.seek(8)

        self.cur_chunk += 1
        self.cpos = 0

    def dump(self):
        self.fileobj.seek(0)
        block_num = 0

        while self.pos < self.fsize:
            print "[BLOCK %s]" % block_num
            # read header
            id1, id2, cm, flg, mtime, xfl, os, xlen = self._read_fields('<BBBBIBBH')
            print 'id1: %s' % id1
            print 'id2: %s' % id2
            print 'cm: %s' % cm
            print 'flg: %s' % flg
            print 'mtime: %s' % mtime
            print 'xfl: %s' % xfl
            print 'os: %s' % os
            print 'xlen: %s' % xlen

            # read subfields
            subpos = 0
            bsize = 0

            while subpos < xlen:
                si1, si2, slen = self._read_fields('<BBH')
                print '    si1: %s' % si1
                print '    si1: %s' % si2
                print '    slen: %s' % slen
                print '    data: [%s]' % slen

                if si1 == 66 and si2 == 67:
                    bsize, = self._read_fields('<H')
                else:
                    self.fileobj.seek(slen, 1)
                    self.pos += slen

                subpos += 6 + slen

            cdata_size = bsize - xlen - 19

            print 'bsize: %s' % bsize
            print 'cdata: [%s]' % (cdata_size)

            self.fileobj.seek(cdata_size, 1)
            self.pos += cdata_size
            crc, isize = self._read_fields('<II')

            print "crc: %s" % crc
            print "isize: %s" % isize
            # print "POS: %s" % self.pos

            block_num += 1

    def _read_fields(self, field_types):
        size = struct.calcsize(field_types)
        self.pos += size
        return struct.unpack(field_types, self.fileobj.read(size))


if __name__ == '__main__':
    print BGZip(sys.argv[1]).dump()
