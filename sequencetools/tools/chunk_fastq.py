#!/usr/bin/env python

import os
import sys
import click
import logging
from signal import signal, SIGPIPE, SIG_DFL
from ..helpers.file_helpers import (return_filehandle, create_directories, 
                                  return_output_handle)
from ..helpers.sequence_helpers import get_seqio_fastq_record

signal(SIGPIPE, SIG_DFL)


def get_chunk(chunks_dir, total_files, gzip_me):
    '''Return new chunk output handle'''
    chunk = '{}/{:06d}.fastq'.format(chunks_dir, total_files)
    if gzip_me:
        chunk += '.gz'
    return return_output_handle(chunk, gzip_me)


def write_chunk(record, chunk, gzip_me):
    '''Formatter for comrpessed and text printing'''
    output = record.format('fastq')
    if gzip_me:
        output = output.encode('utf-8')
    chunk.write(output)


def chunk_fastq(fastq, chunks, chunks_dir, gzip_me):
    '''Chunk FASTQ file.  Output files with chunks reads to chunks_dir
       
       Returns a string with file number and read counts
    '''
    seqio_in = sys.stdin
    fh = ''
    count = 0
    total_reads = 0
    total_files = 1
    create_directories(os.path.abspath(chunks_dir))  # create chunks directory
    chunk = get_chunk(chunks_dir, total_files, gzip_me)
    if not fastq:  # Check STDIN
        for record in get_seqio_fastq_record(seqio_in):  # get SeqIO record
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
        fh = return_filehandle(fastq)
        for record in get_seqio_fastq_record(fh):  # Get SeqIO record
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
    result_str = 'Output {} reads in {} files {} at a time'.format(total_reads,
                                                                   total_files,
                                                                   chunks)
    return result_str


@click.command()
@click.option('--fastq', help='''FASTQ file to chunk, can be compressed''')
@click.option('--chunk_size', help='''Write N reads to file (default:10000)''',
              default=10000)
@click.option('--chunk_dir', 
              help='''Directory to write chunks in (default:./chunks)''', 
              default='./chunks')
@click.option('--gzip_output', is_flag=True,
              help='''Gzip output files''')
@click.option('--log_file', default='./chunk_fastq.log',
              help='''File to write log to.  (default:./chunk_fastq.log)''')
@click.option('--log_level', default='INFO',
    help='''Log level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default:INFO)''')
def main(fastq, chunk_dir, gzip_output, chunk_size, log_file, log_level):
    '''Chunk FASTQ Files.

        cat input*.fastq | chunk_fastq.py

        or

        chunk_fastq.py --fastq input.fastq
    '''
    log_level = getattr(logging, log_level.upper(), logging.INFO)
    msg_format = '%(asctime)s|%(name)s|[%(levelname)s]: %(message)s'
    logging.basicConfig(format=msg_format, datefmt='%m-%d %H:%M',
                        level=log_level)
    log_handler = logging.FileHandler(log_file, mode='w')
    formatter = logging.Formatter(msg_format)
    log_handler.setFormatter(formatter)
    logger = logging.getLogger('chunk_fastq')
    logger.addHandler(log_handler)
    if fastq:
        fastq = os.path.abspath(fastq)
    result = chunk_fastq(fastq, chunk_size, chunk_dir, gzip_output)
    logger.info(result)
        

if __name__ == '__main__':
    main() 
