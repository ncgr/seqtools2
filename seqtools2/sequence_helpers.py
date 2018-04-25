#!/usr/bin/env python

import os
import sys
import select
from Bio import SeqIO


def check_stdin(handle):
    '''Check STDIN using select'''
    if select.select([handle,], [], [], 0.0)[0]:  # use select to check STDIN
        return True  # if True return True
    return False


def get_seqio_record(seq_handle):
    '''Parses a filehandle seq_handle and yields the formatted records
    
       Generator for SeqIO record objects
    '''
    with seq_handle as sopen:
        for record in SeqIO.parse(sopen, 'fasta'):  # iterate with SeqIO
            yield record  # yield each record as it is iterated


if __name__ == '__main__':
    print('Please import!')
    sys.exit(0)
