#!/usr/bin/env python

import os
import sys
import subprocess
import argparse
import logging
import re
from helpers.file_helpers import return_filehandle
from intervaltree import IntervalTree
from time import sleep


parser = argparse.ArgumentParser(description='''

    Detect Discordant Queries for Multiple Reference Sequences

''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--blast_fmt6', metavar = '<blast_fmt6.tbl>', 
required=True,
help='''Blast output to parse.  Format 6 requred.\n\n''')

parser.add_argument('--reference', metavar = '<reference.fasta>',
help='''Required for description incorporation\n\n''')

parser.add_argument('--query_fasta', metavar = '<query.fasta>',
help='''The FASTA file used as the query.  

        Get query length and mark sequences\n\n''')

parser.add_argument('--log_file', metavar = '<FILE>', 
default='./detect_chimeric_alignments.log',
help='''File to write log to.  (default:./detect_chimeric_alignments.log)''')

parser.add_argument('--log_level', metavar = '<LOGLEVEL>', default='INFO',
help='''Log level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default:INFO)''')

parser._optionals.title = "Program Options"
args = parser.parse_args()

log_file = args.log_file

log_level = getattr(logging, args.log_level.upper(), logging.INFO)
msg_format = '%(asctime)s|%(name)s|[%(levelname)s]: %(message)s'
logging.basicConfig(format=msg_format, datefmt='%m-%d %H:%M',
                            level=log_level)
log_handler = logging.FileHandler(
                               log_file, mode='w')
formatter = logging.Formatter(msg_format)
log_handler.setFormatter(formatter)
logger = logging.getLogger('detect_incongruencies')
logger.addHandler(log_handler)


def parse_fasta_headers(reference, ids):
    '''Parse ids and descriptions from FASTA headers of 
    
       reference store in ids dict
    '''
    re_header = re.compile("^>(\S+)\s*(.*)")  # grab id  and desc
    fh = return_filehandle(reference)
    with fh as fopen:
        for line in fopen:
            line = line.rstrip()
            if not line:  # skip blanks
                continue
            if re_header.match(line):  # check for fasta header
                hid = re_header.search(line)
                if hid:
                    logger.debug(hid.groups(0))
                    if isinstance(hid, str):  # check for tuple
                        hid = hid.groups(0)
                    else:
                        rid = hid.groups(0)[0]  # get id portion of header
                        desc = hid.groups(0)[1]
                        ids[rid] = desc
                else:
                    logger.error('Header {} looks odd...'.format(line))
                    sys.exit(1)
                logger.debug(hid)
    return


def check_interval(tree, fields, query_targets, **kwargs):
    '''Check current line from blast fmt6 against existing intervals 
    
       for this reference
    
       if intervals are found to be unique they are added to the tree and the
       
       reference target is assigned, if it is also found to be unique 
       
       for this query
    '''
    if len(fields) < 12:  # Make sure format looks at least rightish
        logger.error('Blast Table should have 12 or more fields')
        sys.exit(1)
    query = fields[0]
    ref = fields[1]
    percent_id = float(fields[2])
    alignment_length = int(fields[3])
    q_start = int(fields[6])
    q_stop = int(fields[7])
    r_start = int(fields[8])
    r_stop = int(fields[9])
    bit_score = float(fields[11])
    sense = '+'
    if q_start > q_stop:  # aligned as reverse complement sense -
        hold = q_start
        q_start = q_stop
        q_stop = hold
        sense = '-'
    query_interval = (q_start, q_stop)
    reference_interval = (r_start, r_stop)
#    data = {'query': query, 'reference': ref, 'q_interval': query_interval,
#            'r_interval': reference_interval, 'sense': sense}
    data = {'q_interval': query_interval, 'r_interval': reference_interval, 
            'sense': sense}
    intervals_envelop = tree.search(q_start, q_stop, strict=True)
    intervals_overlap = tree.search(q_start, q_stop)
    if intervals_envelop:
        logger.debug('Interval {} is enveloped by {}.  Skipping...'.format(
                                                            query_interval,
                                                            intervals_envelop))
        return
    if intervals_overlap:  # Should be option here for HSP combining
        logger.debug('Interval {} overlaps {}.'.format(query_interval,
                                                intervals_overlap))
        for i in intervals_overlap:
            logger.debug('Interval: {},{}'.format(i.begin, i.end))
        return
    tree.addi(q_start, q_stop, data)  # instert into interval tree
    if query not in query_targets:
        query_targets[query] = {}
    if ref not in query_targets[query]:
        query_targets[query][ref] = []
    query_targets[query][ref].append(data)


def detect_chimeras(blast_table, reference, **kwargs):
    '''Main workflow method'''
    blast_table = os.path.abspath(blast_table)
    reference_ids = {}  # store these to attach descriptions later
    query = None  # set None
    query_targets = {}
    count = 0
    out_file = './putative_chimeras.out'
    out_fh = open(out_file, 'w')
    if reference:
        reference = os.path.abspath(reference)
        parse_fasta_headers(reference, reference_ids)
    fh = return_filehandle(blast_table)  # get a filehandle
    with fh as bopen:
        for line in bopen:
            line = line.rstrip()  # remove blank and returns
            fields = line.split('\t')  # split to set values below
            if query != fields[0]:
                if query:
                    if len(query_targets[query]) > 1:
                        logger.info('{}\t{}\t{}'.format(query, 
                                                 query_targets[query], count))
                        out_fh.write('{}\t{}\t{}\n'.format(query, 
                                                  query_targets[query], count))
                tree = IntervalTree()  # init object for interval tree
                count = 0
            query = fields[0]
            ref_id = fields[1]
            if ref_id in reference_ids:  # attach description if true
                fields[1] = ref_id + ' {}'.format(reference_ids[ref_id])
            check_interval(tree, fields, query_targets, **kwargs)
            count += 1
    out_fh.close()
           

if __name__ == '__main__':
    blast_table = args.blast_fmt6
    reference = args.reference
    options = {}  # fill in later for percent identity, length, etc filtering
    detect_chimeras(blast_table, reference, **options)  # main method for workflow
