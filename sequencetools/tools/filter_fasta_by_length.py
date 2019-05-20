#!/usr/bin/env python

import os
import sys
import click
import logging
from signal import signal, SIGPIPE, SIG_DFL
from ..helpers.file_helpers import return_filehandle
from ..helpers.sequence_helpers import get_seqio_fasta_record, check_sequence_length

signal(SIGPIPE, SIG_DFL)


def filter_fasta_by_length(fasta, length, reverse):
    '''Filter FASTA file fasta >= length.

       If reverse, fasta <= length
    '''
    seqio_in = sys.stdin
    fh = ''
    if not fasta:  # Check STDIN
        for record in get_seqio_fasta_record(seqio_in):  # get SeqIO record
            if check_sequence_length(record.seq, length, reverse):  #  length
                print('>{}\n{}'.format(record.description, record.seq))
    else:  # Check FASTA
        fh = return_filehandle(fasta)
        for record in get_seqio_fasta_record(fh):  # Get SeqIO record
            if check_sequence_length(record.seq, length, reverse):  # length
                print('>{}\n{}'.format(record.description, record.seq))


@click.command()
@click.option('--fasta',
    help='''FASTA file to filter, can be compressed''')
@click.option('--length',
    help='''Length Cutoff (default:1000)''', default=1000)
@click.option('--reverse', is_flag=True,
    help='''Filter sequences "<=" instaed of ">="''')
@click.option('--log_file', default='./filter_fasta_by_length.log',
    help='''File to write log to.  (default:./filter_fasta_by_length.log)''')
@click.option('--log_level', default='INFO',
    help='''Log level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default:INFO)''')
def main(fasta, length, reverse, log_file, log_level):
    '''Length Filter for FASTA Files

        cat input.fasta | filter_fasta_by_length.py

        or

        filter_fasta_by_length.py --fasta input.fasta
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
    filter_fasta_by_length(fasta, length, reverse) 


if __name__ == '__main__':
    main()
