#!/usr/bin/env python

import os
import sys
import argparse
import logging
from signal import signal, SIGPIPE, SIG_DFL
from helpers.file_helpers import (return_filehandle, create_directories, 
                                  return_output_handle)
from helpers.sequence_helpers import get_seqio_fasta_record

signal(SIGPIPE, SIG_DFL)

parser = argparse.ArgumentParser(description='''

        Chunk FASTA Files.

        cat input*.fasta | chunk_fasta.py

        or 

        chunk_fasta.py --fasta input.fasta

''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--fasta', metavar = '</path/to/my/fasta.fa>',
help='''FASTA file to chunk, can be compressed''')

parser.add_argument('--chunk_size', metavar = '<INT>', type=int,
help='''Write N reads to file (default:1000)''',
default=1000)

parser.add_argument('--chunk_dir', metavar = '</path/to/chunks>',
help='''Directory to write chunks in (default:./chunks)''',
default='./chunks')

parser.add_argument('--gzip_output', action='store_true',
help='''Gzip output files''')

parser.add_argument('--log_file', metavar = '<FILE>',
default='./chunk_fasta.log',
help='''File to write log to.  (default:./chunk_fasta.log)''')

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
logger = logging.getLogger('chunk_fasta')
logger.addHandler(log_handler)


def get_chunk(chunks_dir, total_files, gzip_me):
    '''Return new chunk output handle'''
    chunk = '{}/{:06d}.fasta'.format(chunks_dir, total_files)
    if gzip_me:
        chunk += '.gz'
    return return_output_handle(chunk, gzip_me)


def write_chunk(record, chunk, gzip_me):
    '''Formatter for comrpessed and text printing'''
    output = record.format('fasta')
    if gzip_me:
        output = output.encode('utf-8')
    chunk.write(output)


def chunk_fasta(fasta, chunks, chunks_dir, gzip_me):
    '''Chunk FASTA file.  Output files with chunks reads to chunks_dir'''
    seqio_in = sys.stdin
    fh = ''
    count = 0
    total_reads = 0
    total_files = 1
    create_directories(os.path.abspath(chunks_dir))  # create chunks directory
    chunk = get_chunk(chunks_dir, total_files, gzip_me)
    if not fasta:  # Check STDIN
        logger.info('Parsing STDIN... Chunking reads in {}s...'.format(chunks))
        for record in get_seqio_fasta_record(seqio_in):  # get SeqIO record
            total_reads += 1
            count += 1
            if count > chunks:  # open new file close old file
                count = 0
                total_files += 1
                chunk.close()
                chunk = get_chunk(chunks_dir, total_files, gzip_me)
                write_chunk(record, chunk, gzip_me)
                continue  # dont write sequence twice
            write_chunk(record, chunk, gzip_me)
    else:  # Check FASTA
        logger.info('Parsing FASTA file... Chunking reads in {}s...'.format(
                                                                      chunks))
        fh = return_filehandle(fasta)
        for record in get_seqio_fasta_record(fh):  # Get SeqIO record
            total_reads += 1
            count += 1
            if count > chunks:  # open new file close old file
                count = 0
                total_files += 1
                chunk.close()
                chunk = get_chunk(chunks_dir, total_files, gzip_me)
                write_chunk(record, chunk, gzip_me)
                continue  # don't write sequence twice
            write_chunk(record, chunk, gzip_me)
    chunk.close()  # close last instance of chunk
    logger.info('Output {} reads in {} files {} at a time'.format(total_reads,
                                                                  total_files,
                                                                  chunks))
            

if __name__ == '__main__':
    fasta = args.fasta
    chunk_dir = os.path.abspath(args.chunk_dir)
    gzip_me = args.gzip_output
    if fasta:
        fasta = os.path.abspath(fasta)
    chunks = args.chunk_size
    chunk_fasta(fasta, chunks, chunk_dir, gzip_me)
