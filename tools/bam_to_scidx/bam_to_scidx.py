"""
bam_to_scidx.py

Input: BAM file
Output: Converted scidx file
"""
import optparse
import os
import shutil
import subprocess
import sys
import tempfile

BUFF_SIZE = 1048576


def get_stderr_exception(tmp_err, tmp_stderr):
    """
    Return a stderr string of reasonable size.
    """
    tmp_stderr.close()
    # Get stderr, allowing for case where it's very large.
    tmp_stderr = open(tmp_err, 'rb')
    stderr_str = ''
    buffsize = BUFF_SIZE
    try:
        while True:
            stderr_str += tmp_stderr.read(buffsize)
            if not stderr_str or len(stderr_str) % buffsize != 0:
                break
    except OverflowError:
        pass
    tmp_stderr.close()
    return stderr_str


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)


parser = optparse.OptionParser()
parser.add_option('-j', '--jar_file', dest='jar_file', type='string', help='BAMtoIDX.jar')
parser.add_option('-b', '--input_bam', dest='input_bam', type='string', help='Input dataset in BAM format')
parser.add_option('-i', '--input_bai', dest='input_bai', type='string', help='Input dataset index')
parser.add_option('-r', '--read', dest='read', type='int', default=0, help='Reads.')
parser.add_option('-o', '--output', dest='output', type='string', help='Output dataset in scidx format')
options, args = parser.parse_args()

tmp_dir = tempfile.mkdtemp(prefix='tmp-scidx-')
tmp_out = tempfile.NamedTemporaryFile().name
tmp_stdout = open(tmp_out, 'wb')
tmp_err = tempfile.NamedTemporaryFile().name
tmp_stderr = open(tmp_err, 'wb')

# Link input BAM file and associated index file into the working directory.
input_data_file_name = "input.bam"
os.link(options.input_bam, os.path.join(tmp_dir, input_data_file_name))
os.link(options.input_bai, os.path.join(tmp_dir, "%s.bai" % input_data_file_name))

cmd = "java -jar %s -b %s -s %d -o %s" % (options.jar_file, input_data_file_name, options.read, options.output)
proc = subprocess.Popen(args=cmd, stdout=tmp_stdout, stderr=tmp_stderr, shell=True, cwd=tmp_dir)
return_code = proc.wait()
if return_code != 0:
    error_message = get_stderr_exception(tmp_err, tmp_stderr)
    stop_err(error_message)

if tmp_dir and os.path.exists(tmp_dir):
    shutil.rmtree(tmp_dir)
