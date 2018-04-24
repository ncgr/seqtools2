#!/usr/bin/env python

import os
import sys
import argparse
import logging
from file_helpers import return_filehandle
from Bio import SeqIO

parser = argparse.ArgumentParser(description='''

        Length Filter for FASTA Files

        cat input.fasta | fasta_filter_by_len.py

        or 

        fasta_filter_by_len.py --fasta input.fasta

''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--fasta', metavar = '</path/to/my/fasta.fa>',
help='''FASTA file to filter, can be compressed''')

parser.add_argument('--length', metavar = '<INT>', type=int,
help='''Length Cutoff (default:1000)''',
default=1000)

parser.add_argument('--reverse', action='store_true',
help='''Filter sequences "<=" instaed of ">="''')

parser.add_argument('--log_file', metavar = '<FILE>',
default='./filter_fasta_by_length.log',
help='''File to write log to.  (default:./filter_fasta_by_length.log)''')

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
logger = logging.getLogger('detect_incongruencies')
logger.addHandler(log_handler)


def check_length(seq, length, reverse):
    '''Accepts a string record and checks to see if it returns true or false

       based on length and reverse
    '''
    if reverse:
        if len(seq) <= length:
            return True
        else:
            return False
    if len(seq) >= length:
        return True
    return False


def get_seqio_record(seq_handle):
    '''Parses a filehandle seq_handle and yields the formatted records'''
    with seq_handle as sopen:
        for record in SeqIO.parse(sopen, 'fasta'):  # iterate with SeqIO
            yield record  # yield each record as it is iterated


def filter_fasta_by_length(fasta, length, reverse):
    '''Filter FASTA file fasta >= length.

       If reverse, fasta <= length
    '''
    seqio_in = sys.stdin
    fh = ''
    if seqio_in:  # Check STDIN
        logger.info('Parsing STDIN...  Checking Sequence Lengths...')
        for record in get_seqio_record(seqio_in):  # get SeqIO record
            if check_length(record.seq, length, reverse):  # check length
                print('>{}\n{}'.format(record.description, record.seq))
    if fasta:  # Check FASTA
        logger.info('Parsing FASTA file...  Checking Sequence Lengths...')
        fh = return_filehandle(fasta)
        for record in get_seqio_record(fh):  # Get SeqIO record
            if check_length(record.seq, length, reverse):  # check length
                print('>{}\n{}'.format(record.description, record.seq))
            

if __name__ == '__main__':
    fasta = args.fasta
    length = args.length
    reverse = args.reverse
    if not fasta and not sys.stdin:
        logger.error('FASTA input file or STDIN required!')
        sys.exit(1)
    if fasta:
        fasta = os.path.abspath(fasta)
    filter_fasta_by_length(fasta, length, reverse)
