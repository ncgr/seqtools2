#!/usr/bin/env python

import os
import sys
import argparse
import logging
from signal import signal, SIGPIPE, SIG_DFL
from helpers.file_helpers import load_targets_file, return_filehandle
from helpers.sequence_helpers import get_seqio_fastq_record, check_sequence_id

signal(SIGPIPE, SIG_DFL)

parser = argparse.ArgumentParser(description='''

        Get a subset of FASTQ sequences from a file by id

        cat input.fastq | get_fastq_by_id.py

        or 

        get_fastq_by_id.py --fastq input.fastq --targets targets.txt

''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--fastq', metavar = '</path/to/my/fastq.fa>',
help='''fastq file to filter, can be compressed''')

parser.add_argument('--targets', metavar='<targets.txt>', required=True,
help='''Targets file, one per line''')

parser.add_argument('--reverse', action='store_true',
help='''Reverses target behavior.  Ignore sequences in targets.txt''')

parser.add_argument('--log_file', metavar = '<FILE>',
default='./get_fastq_by_id.log',
help='''File to write log to.  (default:./get_fastq_by_id.log)''')

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
logger = logging.getLogger('get_fastq_by_id')
logger.addHandler(log_handler)


def get_fastq_by_id(fastq, targets_file, reverse):
    '''Get IDs from targets_file and return fastq records from fastq

       that match the loaded IDs
    '''
    seqio_in = sys.stdin
    fh = ''
    targets = load_targets_file(targets_file)
    if not fastq:  # Check STDIN
        logger.info('Parsing STDIN...  Checking for IDs...')
        for record in get_seqio_fastq_record(seqio_in):  # get SeqIO record
            if check_sequence_id(record.id, targets, reverse):  # check
                sys.stdout.write(record.format('fastq'))
                sys.stdout.flush()
    else:  # Check fastq
        logger.info('Parsing fastq file...  Checking for IDs...')
        fh = return_filehandle(fastq)
        for record in get_seqio_fastq_record(fh):  # Get SeqIO record
            if check_sequence_id(record.id, targets, reverse):  # check
                sys.stdout.write(record.format('fastq'))
                sys.stdout.flush()


if __name__ == '__main__':
    fastq = args.fastq
    targets = args.targets
    reverse = args.reverse
    if fastq:
        fastq = os.path.abspath(fastq)
    if targets:
        targets = os.path.abspath(targets)
    get_fastq_by_id(fastq, targets, reverse)
