#!/usr/bin/env python

import os
import sys
import click
import logging
from signal import signal, SIGPIPE, SIG_DFL
from ..helpers.file_helpers import return_filehandle
from ..helpers.sequence_helpers import get_seqio_fasta_record, check_sequence_length

signal(SIGPIPE, SIG_DFL)


def break_lines(sequence, regions, length):
    '''Breaks lines into strings of length "length" and populate regions list

       this is printed with a join on newline
    '''
    start = 0
    stop = length
    sequence_length = len(sequence)
    while(stop < sequence_length):
        my_region = sequence[start:stop]
        regions.append(my_region)
        start = stop
        stop += length
    my_region = sequence[start:]
    regions.append(my_region)


def format_fasta(fasta, line_length):
    '''Format FASTA file with sequence length line_length.

       will add reheader later
    '''
    fh = sys.stdin
    if fasta:  # Check STDIN
        fh = return_filehandle(fasta)
    for record in get_seqio_fasta_record(fh):  # Get SeqIO record
        regions = []
        break_lines(str(record.seq), regions, line_length)  # build regions for output
        print('>{}\n{}'.format(record.description, '\n'.join(regions)))


@click.command()
@click.option('--fasta',
    help='''FASTA file to filter, can be compressed''')
@click.option('--line_length',
    help='''Length Cutoff (default:80)''', default=80)
@click.option('--log_file', default='./filter_fasta_by_length.log',
    help='''File to write log to.  (default:./filter_fasta_by_length.log)''')
@click.option('--log_level', default='INFO',
    help='''Log level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default:INFO)''')
def main(fasta, line_length, log_file, log_level):
    '''Format FASTA Files

        cat input.fasta | format_fasta.py

        or

        format_fasta.py --fasta input.fasta
    '''
    log_level = getattr(logging, log_level.upper(), logging.INFO)
    msg_format = '%(asctime)s|%(name)s|[%(levelname)s]: %(message)s'
    logging.basicConfig(format=msg_format, datefmt='%m-%d %H:%M',
                        level=log_level)
    log_handler = logging.FileHandler(log_file, mode='w')
    formatter = logging.Formatter(msg_format)
    log_handler.setFormatter(formatter)
    logger = logging.getLogger('filter_fasta_by_length')
    logger.addHandler(log_handler)
    if fasta:  # if not stdin get full path
        fasta = os.path.abspath(fasta)
    format_fasta(fasta, line_length) 


if __name__ == '__main__':
    main()
