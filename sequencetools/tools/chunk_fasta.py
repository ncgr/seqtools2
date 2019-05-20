#!/usr/bin/env python

import os
import sys
import click
import logging
from signal import signal, SIGPIPE, SIG_DFL
from ..helpers.file_helpers import (return_filehandle, create_directories,
                                     return_output_handle)
from ..helpers.sequence_helpers import get_seqio_fasta_record

signal(SIGPIPE, SIG_DFL)


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


def process_filehandle(fh, chunks, chunks_dir, gzip_me, byte_chunks):
    count = 0
    total_reads = 0
    total_files = 0
    create_directories(os.path.abspath(chunks_dir))  # create chunks directory
    chunk = get_chunk(chunks_dir, total_files, gzip_me)
    for record in get_seqio_fasta_record(fh):  # get SeqIO record
        total_reads += 1
        if byte_chunks:  # count is incremented by bytes of sequence
            count += len(record.seq)  # bytes of current sequence
        else:
            count += 1
        if count >= chunks:  # open new file close old file
            count = 0
            chunk.close()
            chunk = get_chunk(chunks_dir, total_files, gzip_me)
            write_chunk(record, chunk, gzip_me)
            total_files += 1
            continue  # dont write sequence twice
        write_chunk(record, chunk, gzip_me)
    chunk.close()  # close last instance of chunk
    result_str = 'Output {} reads in {} files {} at a time'.format(total_reads,
                                                                   total_files,
                                                                   chunks)
    return result_str


def chunk_fasta(fasta, chunks, chunks_dir, gzip_me, byte_chunks):
    '''Chunk FASTA file.  Output files with chunks reads to chunks_dir
    
       if byte_chunks, chunk by bytes.  Will try to put chunks bytes in file.

       Will not split sequences.
    '''
    seqio_in = sys.stdin
    fh = ''
    if not fasta:  # Check STDIN
        return process_filehandle(seqio_in, chunks, chunks_dir, 
                                  gzip_me, byte_chunks)
    else:  # Check FASTA
        fh = return_filehandle(fasta)
        return process_filehandle(fh, chunks, chunks_dir, 
                                  gzip_me, byte_chunks)


@click.command()
@click.option('--fasta', help='''FASTA file to chunk, can be compressed''')
@click.option('--chunk_size', help='''Write N reads to file (default:1000)''',
              default=1000)
@click.option('--chunk_bytes', 
      help='''Try to write N sequence bytes to file.  Keeps sequence intact''')
@click.option('--chunk_dir', 
              help='''Directory to write chunks in (default:./chunks)''',
              default='./chunks')
@click.option('--gzip_output', is_flag=True,
              help='''Gzip output files''')
@click.option('--log_file', default='./chunk_fasta.log',
              help='''File to write log to.  (default:./chunk_fasta.log)''')
@click.option('--log_level', default='INFO',
    help='''Log level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default:INFO)''')
def main(fasta, chunk_dir, chunk_size, gzip_output, 
         chunk_bytes, log_file, log_level):
    '''Chunk FASTA Files.

         cat input*.fasta | chunk_fasta.py

         or

         chunk_fasta.py --fasta input.fasta
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
    byte_chunks = False
    if fasta:
        fasta = os.path.abspath(fasta)
    if chunk_bytes:
        chunk_size = int(chunk_bytes)
        byte_chunks = True
    result = chunk_fasta(fasta, chunk_size, chunk_dir, 
                         gzip_output, byte_chunks)
    logger.info(result)


if __name__ == '__main__':
    main()
