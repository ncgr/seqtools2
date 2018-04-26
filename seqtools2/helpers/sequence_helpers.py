#!/usr/bin/env python

import os
import sys
import select
from Bio import SeqIO


def check_sequence_length(seq, length, reverse):
    '''Accepts a string record and checks to see if it returns true or false
       
       based on length and reverse
    '''
    if reverse:  # make less than
        if len(seq) <= length:
            return True
        else:
            return False
    if len(seq) >= length:  # normal >= check
        return True
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
