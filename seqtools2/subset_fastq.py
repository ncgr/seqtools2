#!/usr/bin/env python

import os
import sys
import argparse
import logging
from signal import signal, SIGPIPE, SIG_DFL
from helpers.file_helpers import return_filehandle
from helpers.sequence_helpers import get_seqio_fastq_record

signal(SIGPIPE, SIG_DFL)

parser = argparse.ArgumentParser(description='''

        Subset FASTQ Files.

        cat input*.fastq | subset_fastq.py

        or 

        subset_fastq.py --fastq input.fastq

''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--fastq', metavar = '</path/to/my/fastq.fq>',
help='''FASTQ file to subset, can be compressed''')

parser.add_argument('--subset', metavar = '<INT>', type=int,
help='''Take every N reads (default:10)''',
default=10)

parser.add_argument('--log_file', metavar = '<FILE>',
default='./subset_fastq.log',
help='''File to write log to.  (default:./subset_fastq.log)''')

parser.add_argument('--log_level', metavar = '<LOGLEVEL>', default='INFO',
help='''Log level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default:INFO)''')

parser._optionals.title = "Program Options"
args = parser.parse_args()

log_file = args.log_file

log_level = getattr(logging, args.log_level.upper(), logging.INFO)
msg_format = '%(asctime)s|%(name)s|[%(levelname)s]: %(message)s'
logging.basicConfig(format=msg_format, datefmt='%m-%d %H:%M',
                    level=log_level)
log_handler = logging.FileHandler(log_file, mode='w')
formatter = logging.Formatter(msg_format)
log_handler.setFormatter(formatter)
logger = logging.getLogger('subset_fastq')
logger.addHandler(log_handler)


def subset_fastq(fastq, subset):
    '''Subset FASTQ file.  Pick 1/subset reads.

       If reverse, fasta <= length
    '''
    seqio_in = sys.stdin
    fh = ''
    count = 0
    total = 0
    if not fastq:  # Check STDIN
        logger.info('Parsing STDIN... Output every {} reads...'.format(subset))
        for record in get_seqio_fastq_record(seqio_in):  # get SeqIO record
            count += 1
            if count == subset:
                count = 0
                total += 1
                print(record.format('fastq'))
    else:  # Check FASTA
        logger.info('Parsing FASTQ file... Output every {} reads...'.format(
                                                                      subset))
        fh = return_filehandle(fastq)
        for record in get_seqio_fastq_record(fh):  # Get SeqIO record
            count += 1
            if count == subset:
                count = 0
                total += 1
                print(record.format('fastq'))
    logger.info('Output {} reads'.format(total))
            

if __name__ == '__main__':
    fastq = args.fastq
    if fastq:
        fastq = os.path.abspath(fastq)
    subset = args.subset
    subset_fastq(fastq, subset)
