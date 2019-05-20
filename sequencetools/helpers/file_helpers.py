#!/usr/bin/env python

import os
import sys
import gzip
import errno
import re
import select


def check_stdin(handle):
    '''Check STDIN using select'''
    if select.select([handle,], [], [], 0.0)[0]:  # use select to check STDIN
        return True  # if True return True
    return False


def create_directories(dirpath):
    '''make directory path'''
    try:
        os.makedirs(dirpath)
    except OSError as e:
        if e.errno != errno.EEXIST:  # ignore if error is exists else raise
            raise


def check_file_type(check_me):
    '''Checks a files type by splitting on "." and getting the last element'''
    suffix = check_me.split('.')[-1].lower()  # get suffix
    fasta_check = re.compile('fa|fasta|fna')
    fastq_check = re.compile('fq|fastq')
    compression_check = re.compile('gz|bgz|xz|zip|bzip')
    if not suffix:
        return False  # no file type
    if compression_check.match(suffix):
        suffix = check_me.split('.')[-2].lower()
    if fastq_check.match(suffix):
        suffix = 'fastq'  # set standard
    elif fasta_check.match(suffix):
        suffix = 'fasta'  # set standard
    return suffix


def return_output_handle(write_me, gzipped):
    '''Returns an open file handle to write.  
    
       If gzipped=True returns gzip handle
    '''
    if gzipped:
        return gzip.open(write_me, 'wb')  # return comrpessed handle
    return open(write_me, 'w')


def return_filehandle(open_me):
    '''get me a filehandle, common compression or text'''
    magic_dict = {
                  b'\x1f\x8b\x08': 'gz'  # only one supported right now
#                  '\x42\x5a\x68': 'bz2',
#                  '\x50\x4b\x03\x04': 'zip'
                 }
    max_bytes = max(len(t) for t in magic_dict)
    with open(open_me, 'rb') as f:  # get read and binary fixes python3 issues
        s = f.read(max_bytes)
    for m in magic_dict:
        if s.startswith(m):
            t = magic_dict[m]  # get type
            if t == 'gz':
                return gzip.open(open_me, 'rt')  # return handle
#            elif t == 'bz2':
#                return bz2.open(open_me)
#            elif t == 'zip':
#                return zipfile.open(open_me)
    return open(open_me)  # return normal handle if not compressed


def load_targets_file(targets_file):
    '''Load targets_file into targets dict and return dict'''
    fh = return_filehandle(targets_file)
    targets = {}  # targets dict to return
    with fh as topen:
        for line in topen:
            line = line.rstrip()
            if not line or line.startswith('#'):  # skip blank and comments
                continue
            targets[line] = 1  # store targets could maybe use a set here
    return targets


if __name__ == '__main__':
    print('Please import to use create_directories and return_filehandle')
    sys.exit(0)
