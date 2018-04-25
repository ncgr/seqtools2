#!/usr/bin/env python

import os
import sys
import argparse
import logging
from time import sleep
from signal import signal, SIGPIPE, SIG_DFL
from sequence_helpers import check_stdin, get_seqio_record
from file_helpers import return_filehandle

signal(SIGPIPE, SIG_DFL) 

parser = argparse.ArgumentParser(description='''

        Get Basic Contiguity, Continuity, and Sequence Content Stats

        cat input.fasta | basic_fasta_stats.py

        or 

        basic_fasta_stats.py --fasta input.fasta

''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--fasta', metavar = '</path/to/my/fasta.fa>',
help='''FASTA file to filter, can be compressed''')

parser.add_argument('--min_gap', metavar = '<INT>',
help="""Minimum length of consecutive N's to consider a gap""")

parser.add_argument('--log_file', metavar = '<FILE>',
default='./basic_fasta_stats.log',
help='''File to write log to.  (default:./basic_fasta_stats.log)''')

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
logger = logging.getLogger('basic_fasta_stats')
logger.addHandler(log_handler)


def get_N50(contigs, total):
    c = 0
    n = 0
    for l in reversed(contigs):
        c += int(l)
        n += 1
        if c >= float(total)/2:
            return l
    return 0


def get_gc(gc, total):
    if not total or not gc:
        return 0
    pgc = float(gc)/int(total)
    return pgc


def compile_metrics(metrics, lengths, bases):
    lengths['contigs'] = sorted(lengths['contigs'])
    lengths['gaps'] = sorted(lengths['gaps'])
    lengths['scaffolds'] = sorted(lengths['scaffolds'])
    lengths['total'] = sorted(lengths['total'])
    if lengths['contigs']:
        metrics['maxcontig'] = lengths['contigs'][-1]
        metrics['mincontig'] = lengths['contigs'][0]
        if metrics['contigbases']:
            metrics['contigN50'] = get_N50(lengths['contigs'], 
                                           metrics['contigbases'])
    if lengths['gaps']:
        metrics['maxgap'] = lengths['gaps'][-1]
        metrics['mingap'] = lengths['gaps'][0]
        if metrics['gapbases']:
            metrics['gapN50'] = get_N50(lengths['gaps'], metrics['gapbases'])
    if lengths['scaffolds']:
        metrics['maxscaffold'] = lengths['scaffolds'][-1]
        metrics['minscaffold'] = lengths['scaffolds'][0]
        if metrics['scaffoldbases']:
            metrics['scaffoldN50'] = get_N50(lengths['scaffolds'], 
                                             metrics['scaffoldbases'])
    metrics['allbases'] = bases['total']
    metrics['N50'] = get_N50(lengths['total'], bases['total'])


def basic_fasta_stats(fasta):
    if not fasta:  # Assume STDIN
        fasta = sys.stdin
    else:
        fasta = return_filehandle(fasta)
    bases = {'A' : 0, 'a' : 0, 'C' : 0, 'c' : 0,
             'T' : 0, 't' : 0, 'G' : 0, 'g' : 0,
             'N' : 0, 'n' : 0, 'IUPAC' : 0, 'total' : 0}
    metrics = {'N50' : 0, 'N90' : 0, 'L50' : 0,
               'contigN50' : 0, 'scaffoldN50' : 0,
               'gaps' : 0, 'gapN50' : 0, 'maxgap' : 0, 'mingap' : 0,
               'contigs' : 0, 'scaffolds' : 0, 'records' : 0,
               'maxscaffold' : 0, 'minscaffold' : 0,
               'maxcontig' : 0, 'mincontig' : 0,
               'contigbases' : 0, 'scaffoldbases' : 0, 'gapbases' : 0,
               'allbases' : 0, 'pgc' : 0}
    min_gap = args.min_gap
    lengths = {'scaffolds' : [], 'gaps' : [], 'contigs' : [], 'total' : []}
    scheck = 0
    gcheck = 0
    ccheck = 0
    glen = 0
    clen = 0
    slen = 0
    length = 0
    i = 0
    for record in get_seqio_record(fasta):
        metrics['records'] += 1
        seq = record.seq
        if seq:
            seq = seq.upper()
            length = len(seq)
            bases['total'] += length
            lengths['total'].append(length)
            while True:
                i = seq.find('N', i)
                if i < 0:
                    if not scheck:
                        lengths['contigs'].append(length)
                        metrics['contigs'] += 1
                        metrics['contigbases'] += length
                        for b in seq:
                            if b in bases:
                                bases[b] += 1
                            else:
                                bases['IUPAC'] += 1
                    else:
                        clen = length - ccheck
                        lengths['contigs'].append(clen)
                        metrics['contigs'] += 1
                        metrics['contigbases'] += clen
                        for b in seq[ccheck:]:
                            if b in bases:
                                bases[b] += 1
                            else:
                                bases['IUPAC'] += 1
                    break
                gcheck = i
                while i < length and seq[i] == 'N':
                    i += 1
                glen = i - gcheck
                if glen >= min_gap:
                    scheck = 1
                    clen = gcheck - ccheck
                    lengths['gaps'].append(glen)
                    bases['N'] += glen
                    metrics['gaps'] += 1
                    metrics['gapbases'] += glen
                    lengths['contigs'].append(clen)
                    metrics['contigs'] += 1
                    metrics['contigbases'] += clen
                    for b in seq[ccheck:gcheck]:
                        if b in bases:
                            bases[b] += 1
                        else:
                            bases['IUPAC'] += 1
                    ccheck = i
            if scheck:
                metrics['scaffolds'] += 1
                lengths['scaffolds'].append(length)
                metrics['scaffoldbases'] += length
        scheck = 0
        gcheck = 0
        ccheck = 0
        glen = 0
        clen = 0
        i = 0
        pos = 0
    compile_metrics(metrics, lengths, bases)
    return metrics


if __name__ == '__main__':
    fasta = args.fasta
    if fasta and check_stdin(sys.stdin):
        logger.error('Can only provide STDIN or FASTA file, not both!')
        sys.exit(1)
    if fasta:
        fasta = os.path.abspath(fasta)
    print(basic_fasta_stats(fasta))