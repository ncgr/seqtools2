#!/usr/bin/env python

import os
import sys
import gzip
import errno

def create_directories(dirpath):
    '''make directory path'''
    try:
        os.makedirs(dirpath)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def return_filehandle(open_me):
    '''get me a filehandle, common compression or text'''
    magic_dict = {
                  b'\x1f\x8b\x08': 'gz'
#                  '\x42\x5a\x68': 'bz2',
#                  '\x50\x4b\x03\x04': 'zip'
                 }
    max_bytes = max(len(t) for t in magic_dict)
    with open(open_me, 'rb') as f:
        s = f.read(max_bytes)
    for m in magic_dict:
        if s.startswith(m):
            t = magic_dict[m]
            if t == 'gz':
                return gzip.open(open_me, 'rt')
#            elif t == 'bz2':
#                return bz2.open(open_me)
#            elif t == 'zip':
#                return zipfile.open(open_me)
    return open(open_me)


if __name__ == '__main__':
    print('Please import to use create_directories and return_filehandle')
    sys.exit(0)
