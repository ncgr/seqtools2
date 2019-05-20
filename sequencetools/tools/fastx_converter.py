#!/usr/bin/env python

import os
import sys
import re
import click
import logging
from signal import signal, SIGPIPE, SIG_DFL
from ..helpers.file_helpers import return_filehandle, check_file_type
from ..helpers.sequence_helpers import get_seqio_fastx_record


signal(SIGPIPE, SIG_DFL)


def record_to_stdout(record, output_type, quality):
    '''Takes a record, the output_type and the quality and writes to stdout.

       The quality is assigned for fasta to fastq.
    '''
    if output_type == 'fasta':
        sys.stdout.write(record.format(output_type))  # can write normally
        sys.stdout.flush()
    elif output_type == 'fastq':
        length = len(record.seq)
        qualities = [quality] * length
        record.letter_annotations["solexa_quality"] = qualities
        sys.stdout.write(record.format(output_type))  # can write normally
        sys.stdout.flush()


def fastx_converter(input_file, input_type, output_type, quality):
    '''Convert input_file or stdin fasta to fastq or fastq to fasta 
    
       based on input_type
    '''
    seqio_in = sys.stdin
    fh = ''
    if not input_file:  # Check STDIN
        for record in get_seqio_fastx_record(seqio_in, input_type):  # Generatea
            record_to_stdout(record, output_type, quality)
    else:  # Check file
        input_file = os.path.abspath(input_file)
        fh = return_filehandle(input_file)
        for record in get_seqio_fastx_record(fh, input_type):  # Generator
            record_to_stdout(record, output_type, quality)


@click.command()
@click.option('--input_file',
              help='''Input file, fasta or fastq, can be compressed''')
@click.option('--input_type', required=True,
              help='''Input file type.  fasta or fastq''')
@click.option('--output_quality', default=40,
   help='''Quality to assign if converting from fasta to fastq (default:40)''')
@click.option('--log_file', default='./fastx_converter.log',
             help='''File to write log to.  (default:./fastx_converter.log)''')
@click.option('--log_level', default='INFO',
    help='''Log level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default:INFO)''')
def main(input_file, input_type, output_quality, log_file, log_level):
    '''Convert FASTA to FASTQ or FASTQ to FASTA

        cat input.[fa|fq] | fastx_converter.py --input_type <fasta/fastq>

        or

        fastx_converter.py --input input.[fa|fq] --input_type <fasta/fastq>
    '''
    log_level = getattr(logging, log_level.upper(), logging.INFO)
    msg_format = '%(asctime)s|%(name)s|[%(levelname)s]: %(message)s'
    logging.basicConfig(format=msg_format, datefmt='%m-%d %H:%M',
                        level=log_level)
    log_handler = logging.FileHandler(log_file, mode='w')
    formatter = logging.Formatter(msg_format)
    log_handler.setFormatter(formatter)
    logger = logging.getLogger('fastx_converter')
    logger.addHandler(log_handler)
    fasta_check = re.compile('fa|fasta|fna')
    fastq_check = re.compile('fq|fastq')
    output_type = ''
    if fastq_check.match(input_type):  # check input type set fastq
        input_type = 'fastq'
        output_type = 'fasta'
    elif fasta_check.match(input_type):  # check input type set fasta
        input_type = 'fasta'
        output_type = 'fastq'
    else:  # exit if not fasta or fastq
        logger.error('Input type: {} cannot be processed'.format(input_type))
        sys.exit(1)
    if input_file:
        input_type_check = check_file_type(input_file)
        if input_type_check != input_type:  # file doesnt look like input_type
            logger.error('Type mismatch, input_type:{}, file:{}'.format(
                                                             input_type,
                                                             input_type_check))
            sys.exit(1)
    fastx_converter(input_file, input_type, output_type, output_quality)


if __name__ == '__main__':
    main()
