#!/usr/bin/env python

import os
import sys
import click
import logging
from signal import signal, SIGPIPE, SIG_DFL
from ..helpers.file_helpers import return_filehandle
from ..helpers.sequence_helpers import get_seqio_fastq_record

signal(SIGPIPE, SIG_DFL)


def subset_fastq(fastq, subset):
    '''Subset FASTQ file.  Pick 1/subset reads.

       If reverse, fasta <= length
    '''
    seqio_in = sys.stdin
    fh = ''
    count = 0
    total = 0
    if not fastq:  # Check STDIN
        for record in get_seqio_fastq_record(seqio_in):  # get SeqIO record
            count += 1
            if count == subset:
                count = 0
                total += 1
                sys.stdout.write(record.format('fastq'))
                sys.stdout.flush()
    else:  # Check FASTA
        fh = return_filehandle(fastq)
        for record in get_seqio_fastq_record(fh):  # Get SeqIO record
            count += 1
            if count == subset:
                count = 0
                total += 1
                sys.stdout.write(record.format('fastq'))
                sys.stdout.flush()
    return 'Output {} reads'.format(total)


@click.command()            
@click.option('--fastq',
              help='''FASTQ file to subset, can be compressed''')
@click.option('--subset', metavar = '<INT>',
              help='''Take every N reads (default:10)''', default=10)
@click.option('--log_file', metavar = '<FILE>', default='./subset_fastq.log',
              help='''File to write log to.  (default:./subset_fastq.log)''')
@click.option('--log_level', default='INFO',
    help='''Log level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default:INFO)''')
def main(fastq, subset, log_file, log_level):
    '''Subset FASTQ Files.

        cat input*.fastq | subset_fastq.py

        or

        subset_fastq.py --fastq input.fastq
    '''
    log_level = getattr(logging, log_level.upper(), logging.INFO)
    msg_format = '%(asctime)s|%(name)s|[%(levelname)s]: %(message)s'
    logging.basicConfig(format=msg_format, datefmt='%m-%d %H:%M',
                        level=log_level)
    log_handler = logging.FileHandler(log_file, mode='w')
    formatter = logging.Formatter(msg_format)
    log_handler.setFormatter(formatter)
    logger = logging.getLogger('subset_fastq')
    logger.addHandler(log_handler)
    if fastq:
        fastq = os.path.abspath(fastq)
    logger.info(subset_fastq(fastq, subset))


if __name__ == '__main__':
    main()
