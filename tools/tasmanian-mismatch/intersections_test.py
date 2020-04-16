#!/bin/python
'''
    input sam file from stdin
    bedfile as one of the arguments (chromosome, start, end - delimiter=tab)

'''
import os, sys, gzip, re
import numpy as np
#from pathlib import Path
#sys.path.append(str(Path(__file__).parent.parent.parent) + '/utils/')
#sys.path.append(sys.path[0] + '/utils')
try:
    from tasmanian.utils.sam_reads import reads
    from tasmanian.utils.utils import *

except Exception as e:
    # Either tests or base_dir, it's downstream of ../tasmanian/tasmanian/
    p = os.path.abspath(os.path.dirname(__file__))
    #p = re.search("(.*tasmanian/tasmanian/).*",p).group(1)
    p_start = [i for i in re.finditer('/tasmanian',p)][-1].end()
    p = p[:p_start]
    utils_path = p + '/utils'
    sys.path = [utils_path] + sys.path

    #p = re.search("(.*tasmanian/tasmanian/).*",p).group(1)
    #utils_path = p + 'utils'
    #sys.path = [utils_path] + sys.path
    from sam_reads import reads
    from utils import *



def main():

    HELP = '''
    \t\tsamtools view <bam_file> | python -b <bed_file/bedGraph> -o <output.table>\n
    \t\tthe bedGraph file should contain 3 or more columns separated by tabs:
    \t\t------------------------------------------
    \t\tchrI    850     879     +       L1P5    LINE    L1

    '''

    if len(sys.argv)==1: exit('\n-h|--help\n') # there should be at least one argument = '--reference-genome'


    # initialize lists to contain reads (sam) and statistics (e.g. length of intersections)
    statistics = []
    sam_output = []
    sam_output = []

    # load global arguments
    out_prefix = ''
    debug = False

    for n,i in enumerate(sys.argv):
        if i in ['--bed','-b','--bed-file']:
            bedfile = sys.argv[n+1]
        if i in ['--output', '-o']:
            out_prefix = sys.argv[n+1] + "."
        if i in ['-h','--help']:
            print(HELP)
            sys.exit(1)
        if i in ['-d','--debug']:
            debug = True

    if debug:
        # define the name of the logging file for debugging
        logFileName = out_prefix + 'intersections.log'

        # if debugging create this logfile    
        logging.basicConfig(filename = logFileName,
                            format = '%(asctime)s %(message)s',
                            filemode = 'w')
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)

    # load bed_file
    try:
        bed, bed_lens, total_bed_lens, bed_other_info = read_bed(bedfile)
        if debug:
            logger.info('bedfile {} was succesfully read'.format(bedfile))
    except Exception as e:
        if debug:
            logger.error('{} happened in excecution of read_bed in main()'.format(str(e)))
        exit("there was a problem reading {}. Make sure is tab delimited and all columns and rows are correct.")

    # to avoid going over the entire chromosome on each read,
    # keep an updated index of the bedfile to start from
    total_beds = 0 
    last_chrom = 'chr0'
    finito=False # know when we are finished with the bed file

    # define a table of flags of proper paired reads
    proper_flags = {
        99: 'first fwd',
        147:'second fwd',
        83: 'first rev',
        163:'second rev'
    }

    # buffer to keep reads until the pair is found, so the paired reads can be analyzed together.
    buffer = {}
    n_tester = 0

    # read bam file from stdin and 
    for line in sys.stdin:

        # parse the header UNmodified as this script could be used to save the intersections file.
        # without being parsed into tasmanian
        if line[0]=='@': 
            sys.stdout.write(line)
            continue

        line = line.strip('\n')
            
        # In case this is a regular read
        n_tester +=1

        if len(line) < 50:  # avoid potential empty lines at the end.
            if debug:
                logger.warning('line {} had less than 50 characters'.format(n_tester))
            continue 

        # instantiate object "current_read"
        try:
            line = line.split('\t')
            _id, flag, chrom, start, mapq, cigar, _2, _3, tlen, seq, phred = line[:11]
            tags = '\t'.join(line[11:])  # This will be added to the otuput sam file

            # check that bam file is sorted. Use first 200 reads for that
            # current_read.start is the previous coordinate. Should be lower than the start of the 
            # current read (to be assigned).
            if n_tester>1 and n_tester<500:
                if current_read.chrom == chrom and flag in ['99','163'] and current_read.start > int(start):
                    sys.stdout.write( '\x1b[5;37;41m' + 'Error:' + '\x1b[0m' + 'BAM file should be sorted!!\n')
                    sys.exit(1)

            current_read = reads(_id, flag, chrom, start, mapq, cigar, _2, _3, tlen, seq, phred, tags)

            #print(current_read.__dict__)        

            # assume read is not paired yet
            paired_read = None # assume there is no paired_read yet

        except Exception as e:
            if debug:
                logger.error('read {} could\'t be loaded properly'.format(n_tester))
            continue            

        # only consider uniquely mapped and proper pair
        if current_read.flag not in proper_flags:
            continue

        # if chromosome not in bed, read is non overlapping
        if current_read.chrom not in bed:
            sam_output.append(current_read.print('original'))
            continue

        # new chromosome -> zero bed_index!
        if current_read.chrom != last_chrom:
            if debug:
                #logger.critical('This should ONLY happen ONCE!!! {} != {}'.format(current_read.chrom, last_chrom))
                logger.info('still have {} reads in buffer for chromosome {}. Thrown to trash  --  {}'.format(len(buffer), last_chrom, n_tester))
            skip_chrom=False                    # in case this was on for the previous chromosome
            bed_index=0                         # new chromosome, new bed_index
            buffer = {}                         # restart buffer for the new chromosome
            if total_beds >= total_bed_lens:    # Have we finished reading the bed file?
                finito = True
            last_chrom = current_read.chrom


        # Check-point ====================================================================================
        '''
            If read is first on the pair, save it to memory (buffer) until the paired read is found.
            If read is second on the pair, we already have the first in the memory, in a list of first
            reads (sometimes they are 100s of lines appart) and we compute all what follows for both reads.
            99: first -> I prefer not to rely on this. Exceptions could be many more than I think in 
            different samples. 
        '''
        if current_read._id in buffer:
            paired_read = buffer.pop(current_read._id) #, "None")   error is better than None here
        
        # sanity check
        if paired_read != None: # TESTED
            if current_read.chrom != paired_read.chrom:
                if debug:
                    logger.warning('current_read and paired_read have different chromosomes = chimeras')
                continue

        # Bam Block ====================================================================================== 
        ''' In THIS block bam and bed are updated to same genomic region if there are gaps between the two. '''

        # if finished looping over the bed file or bed is ahead of bam
        if finito or skip_chrom or current_read.end <= bed[chrom][bed_index,0]:
            current_read.category = 1


        # if read is ahead of bed fragment, advance index
        else: # don't get in here if no more bed fragments for that chrom or bed is finished
            while current_read.start >= bed[chrom][bed_index,1]:

                if bed_index == bed_lens[chrom]-1: # If here, bed file finished for this chrome
                    skip_chrom = True
                    break
                else:
                    bed_index +=1
                    total_beds+=1


        # If: 1.bam was lower than bed  2.we skipped bed with no bam coverage 3.we got to a bed that's after 
        # the bam region 4. Even though current_read.category is None, no "ab" (see the following section) 
        # will be assigned. Hence, we can already assign a category=1 here.
        if current_read.end <= bed[chrom][bed_index,0] or skip_chrom:
            current_read.category=1
        
        # We can already assign bed_id       
        current_read.bed_id = bed_index
        current_read.bed_extra_info = bed_other_info[chrom][bed_index]


        # Block for intersecting reads ==================================================================== 
        if current_read.category == None: # I have not assigned 1 to it! -> it intercepts a bed region!

            # capture the intersection into a flag called 
            # 1 read      ------------------- 
            #   bed        a  =========  b
            # 2 read              -----------
            # 3 read      --------
            # 4 read           ------
            ''' 
            cases 2 and 3 with only 1 base in common (...a>=0) is equivalent
            to no intersection since I already added a number to the bed-end
            to include the upper value and the read goes from zero and should
            not include the upper value. a=0 and b=0 is included in case 4.
            '''
            a0 = bed[chrom][bed_index][0] - current_read.start
            b0 = current_read.end - bed[chrom][bed_index][1]
            a = a0 if a0>=0 else 0
            b = -b0 if b0>0 else current_read.seq_len
            
            # add feature to read to classify it later
            current_read.category_positions = [a,b] 
            
            # mask intersections with lower-case letters instead of "N" --> I end up with ~half # reads
            #intersect_size = len(current_read.seq[a:b])
            #current_read.masked_seq = current_read.seq[:a] + ''.join(['N'] * intersect_size) + current_read.seq[b:]
            #current_read.intersect_seq = ''.join(['N'] * a) + current_read.seq[a:b] + ''.join(['N'] * b0)
    
            
            intersect = current_read.seq[a:b]
            current_read.masked_seq = current_read.seq[:a] + intersect.lower() + current_read.seq[b:]
            current_read.intersect_seq = current_read.seq[:a].lower() + intersect + current_read.seq[b:].lower()
            b = b if b>=0 else current_read.seq_len+b
            current_read.junction = '{}.{};{}.{}'.format(a,b,bed[chrom][bed_index][0], bed[chrom][bed_index][1])
            current_read.complement = current_read.seq_len - b + a   # len - (b-a) 

            #print(current_read.masked_seq)


            # cigar is to correlate clips with intersections. It might happen that most 
            # intersections are not correctly mapped and therefore, are softcliped
            current_read.expand_cigar()     # previously called CIGAR
            intersect_size=None

        # Current_read.category is already assigned ========================================================
        # read_category will be included as a flag at the end of each read in bam file.
        # we nead both paired reads to calculate this flag.
        if paired_read == None: 
            buffer[current_read._id] = current_read

        else:  
            current_read.category = assign_category(current_read, paired_read)
            paired_read.category = current_read.category # category is shared by the paired reads

            # Here include how many bases in the complementary region for both reads (this will also give us 
            # how many in the intersection regions indirectly)

            current_read.complement = current_read.complement + paired_read.complement
            paired_read.complement = current_read.complement
        
            sam_output.append(paired_read.print('masked'))   # if masked_seq==None goes for 'original' by default.
            sam_output.append(current_read.print('masked'))  # same here!
        

            #for READ in [paired_read, current_read]:
            #    if READ.intersect_seq != None:
            #        if READ.masked_seq == None: logger.critical('This should not be happening, masked is None and intersect is not...?')
            #        sam_output.append(READ.print('intersect'))


    # save files:
    #for mtx, fle in zip([sam_output, sam_output], ['masked','intersections']):
    #    with gzip.open(out_prefix + fle + '.sam.gz', 'wb') as f:
    #        f.write('\n'.join(mtx).encode())
    # 
    #    print(fle, ' written!')
    #    #del mtx
    #    #del f
    
    '''
    # pipe into 2 subprocess that will run on parallel 
    # from https://gist.github.com/waylan/2353749
    from subprocess import Popen, PIPE
    import errno
    
    p = Popen(['run_tasmanian','-r','/mnt/galaxy/data/genome/grch38_full/bwa/GRCh38_full_analysis_set_plus_decoy_hla.fa'], stdin=PIPE)
    to_print = '\n'.join(sam_output)

    try:
        p.stdin.write(b'to_print')
    except IOError as e:
        if e.errno == errno.EPIPE or e.errno == errno.EINVAL:
            exit('Oh OH!!!')
        else:
            raise

    p.stdin.close()
    p.wait()
    '''
    

    import socket, errno
    
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)

    #try:
    #    sys.stdout.write('\n'.join(sam_output))
    #except socket.error as e:
    #    print('A socket error')
    #except IOError as e:
    #    if e.errno == errno.EPIPE:
    #        print('EPIPE error')
    #    else:
    #        print('something else...')
    sys.stdout.write('\n'.join(sam_output))



if __name__=='__main__': main()
