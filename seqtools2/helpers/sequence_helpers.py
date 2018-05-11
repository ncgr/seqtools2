#!/usr/bin/env python

import os
import sys
import select
from Bio import SeqIO


def check_sequence_id(seq_id, targets, reverse):
    '''Checks seq_id string to see if it is in targets.

       if reverse check to see that it is NOT in targets.
    '''
    if reverse:  # not in targets
        if seq_id not in targets:
            return True
        else:
            return False
    if seq_id in targets:  # in targets
        return True
    return False


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


def get_seqio_fasta_record(seq_handle):
    '''Parses a fasta filehandle seq_handle and yields the formatted records
    
       Generator for SeqIO record objects
    '''
    with seq_handle as sopen:
        for record in SeqIO.parse(sopen, 'fasta'):  # iterate with SeqIO
            yield record  # yield each record as it is iterated


def get_seqio_fastq_record(seq_handle):
    '''Parses a fasta filehandle seq_handle and yields the formatted records
    
       Generator for SeqIO record objects
    '''
    with seq_handle as sopen:
        for record in SeqIO.parse(sopen, 'fastq'):  # iterate with SeqIO
            yield record  # yield each record as it is iterated


def get_seqio_fastx_record(seq_handle, file_type):
    '''Takes file type and the sequence filehandle
    
       Generator for general SeqIO records lets Bio handle exceptions
    '''
    with seq_handle as sopen:
        for record in SeqIO.parse(sopen, file_type):
            yield record  # yeild SerIO record object.  Its a set.


if __name__ == '__main__':
    print('Please import!')
    sys.exit(0)
