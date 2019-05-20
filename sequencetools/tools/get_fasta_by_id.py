#!/usr/bin/env python

import os
import sys
import click
import logging
from signal import signal, SIGPIPE, SIG_DFL
from ..helpers.file_helpers import load_targets_file, return_filehandle
from ..helpers.sequence_helpers import get_seqio_fasta_record, check_sequence_id

signal(SIGPIPE, SIG_DFL)


def get_fasta_by_id(fasta, targets_file, reverse):
    '''Get IDs from targets_file and return FASTA records from fasta

       that match the loaded IDs
    '''
    seqio_in = sys.stdin
    fh = ''
    targets = load_targets_file(targets_file)
    if not fasta:  # Check STDIN
        for record in get_seqio_fasta_record(seqio_in):  # get SeqIO record
            if check_sequence_id(record.id, targets, reverse):  # check
                print('>{}\n{}'.format(record.description, record.seq))
    else:  # Check FASTA
        fh = return_filehandle(fasta)
        for record in get_seqio_fasta_record(fh):  # Get SeqIO record
            if check_sequence_id(record.id, targets, reverse):  # check
                print('>{}\n{}'.format(record.description, record.seq))


@click.command()
@click.option('--fasta', help='''FASTA file to filter, can be compressed''')
@click.option('--targets', required=True, 
              help='''Targets file, one per line''')
@click.option('--reverse', is_flag=True,
         help='''Reverses target behavior.  Ignore sequences in targets.txt''')
@click.option('--log_file', default='./get_fasta_by_id.log',
         help='''File to write log to.  (default:./get_fasta_by_id.log)''')
@click.option('--log_level', default='INFO',
    help='''Log level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default:INFO)''')
def main(fasta, targets, reverse, log_file, log_level):
    '''Get a subset of FASTA sequences from a file by id

        cat input.fasta | get_fasta_by_id.py --targets targets.txt

        or

        get_fasta_by_id.py --fasta input.fasta --targets targets.txt
    '''
    log_level = getattr(logging, log_level.upper(), logging.INFO)
    msg_format = '%(asctime)s|%(name)s|[%(levelname)s]: %(message)s'
    logging.basicConfig(format=msg_format, datefmt='%m-%d %H:%M',
                        level=log_level)
    log_handler = logging.FileHandler(log_file, mode='w')
    formatter = logging.Formatter(msg_format)
    log_handler.setFormatter(formatter)
    logger = logging.getLogger('get_fasta_by_id')
    logger.addHandler(log_handler)
    if fasta:  # get full path to input file
        fasta = os.path.abspath(fasta)
    if targets:  # get full path to targets
        targets = os.path.abspath(targets)
    get_fasta_by_id(fasta, targets, reverse)


if __name__ == '__main__':
    main()
