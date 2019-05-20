#!/usr/bin/env python

import os
import sys
import click
import logging
from signal import signal, SIGPIPE, SIG_DFL
from ..helpers.file_helpers import load_targets_file, return_filehandle
from ..helpers.sequence_helpers import get_seqio_fastq_record, check_sequence_id

signal(SIGPIPE, SIG_DFL)


def print_record(record):
    '''Formatter for comrpessed and text printing'''
    output = record.format('fastq')
    print(output)


def get_fastq_by_id(fastq, targets_file, reverse):
    '''Get IDs from targets_file and return FASTQ records from fastq

       that match the loaded IDs
    '''
    seqio_in = sys.stdin
    fh = ''
    targets = load_targets_file(targets_file)
    if not fastq:  # Check STDIN
        for record in get_seqio_fastq_record(seqio_in):  # get SeqIO record
            if check_sequence_id(record.id, targets, reverse):  # check
                sys.stdout.write(record)
                sys.stdout.flush()
    else:  # Check FASTQ
        fh = return_filehandle(fastq)
        for record in get_seqio_fastq_record(fh):  # Get SeqIO record
            if check_sequence_id(record.id, targets, reverse):  # check
                sys.stdout.write(record)
                sys.stdout.flush()


@click.command()
@click.option('--fastq', help='''FASTQ file to filter, can be compressed''')
@click.option('--targets', required=True, 
              help='''Targets file, one per line''')
@click.option('--reverse', is_flag=True,
         help='''Reverses target behavior.  Ignore sequences in targets.txt''')
@click.option('--log_file', default='./get_fastq_by_id.log',
         help='''File to write log to.  (default:./get_fastq_by_id.log)''')
@click.option('--log_level', default='INFO',
    help='''Log level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default:INFO)''')
def main(fastq, targets, reverse, log_file, log_level):
    '''Get a subset of FASTQ sequences from a file by id

        cat input.fastq | get_fastq_by_id.py --targets targets.txt

        or

        get_fastq_by_id.py --fastq input.fastq --targets targets.txt
    '''
    log_level = getattr(logging, log_level.upper(), logging.INFO)
    msg_format = '%(asctime)s|%(name)s|[%(levelname)s]: %(message)s'
    logging.basicConfig(format=msg_format, datefmt='%m-%d %H:%M',
                        level=log_level)
    log_handler = logging.FileHandler(log_file, mode='w')
    formatter = logging.Formatter(msg_format)
    log_handler.setFormatter(formatter)
    logger = logging.getLogger('get_fastq_by_id')
    logger.addHandler(log_handler)
    if fastq:  # get full path to input file
        fastq = os.path.abspath(fastq)
    if targets:  # get full path to targets
        targets = os.path.abspath(targets)
    get_fastq_by_id(fastq, targets, reverse)


if __name__ == '__main__':
    main()
