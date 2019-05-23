

import os
import sys
import click
import json
import logging
from collections import OrderedDict
from time import sleep
from signal import signal, SIGPIPE, SIG_DFL
from ..helpers.sequence_helpers import get_seqio_fasta_record
from ..helpers.file_helpers import return_filehandle, check_stdin

signal(SIGPIPE, SIG_DFL) 


def get_N50(lengths, total):
    '''Calculates the N50 of the list lengths using a running sum and
    
       the provided total sum of lengths total
    '''
    total_length = 0  # will sum bases for all 
    number = 0
    for l in reversed(lengths):
        total_length += int(l)
        number += 1
        if total_length >= float(total)/2:  # N50 is this length l
            return l
    return 0  # if empty list


def get_N90(lengths, total):
    '''Calculates the N50 of the list lengths using a running sum and
    
       the provided total sum of lengths total
    '''
    total_length = 0  # will sum bases for all 
    number = 0
    for l in reversed(lengths):
        total_length += int(l)
        number += 1
        if total_length >= float(total)*0.9:  # N90 is this length l
            return l
    return 0  # if empty list


def get_mean(lengths):
    '''Calculates the N50 of the list lengths using a running sum and
    
       the provided total sum of lengths total
    '''
    total_length = 0  # will sum bases for all 
    number = 0
    mean = 0
    for l in lengths:
        total_length += int(l)
        number += 1
    if number:
        mean = total_length/number
        return round(mean)
    return 0  # if empty list


def get_gc(gc, total):
    '''Calculates the percentage gc'''
    if not total or not gc:
        return 0
    pgc = float(gc)/int(total)  # get percent GC
    return round(pgc)


def compile_metrics(metrics, lengths, bases):
    '''Fill the metrics dictionary with the results from lengths and bases'''
    lengths['contigs'] = sorted(lengths['contigs'])
    lengths['gaps'] = sorted(lengths['gaps'])
    lengths['scaffolds'] = sorted(lengths['scaffolds'])
    lengths['total'] = sorted(lengths['total'])
    if lengths['contigs']:
        metrics['maxcontig'] = lengths['contigs'][-1]
        metrics['mincontig'] = lengths['contigs'][0]
        if metrics['contigbases']:  # contig N50
            metrics['contigN50'] = get_N50(lengths['contigs'], 
                                           metrics['contigbases'])
            metrics['contigN90'] = get_N90(lengths['contigs'], 
                                           metrics['contigbases'])
            metrics['meancontig'] = get_mean(lengths['contigs'])
    if lengths['gaps']:
        metrics['maxgap'] = lengths['gaps'][-1]
        metrics['mingap'] = lengths['gaps'][0]
        if metrics['gapbases']:  # gap N50
            metrics['gapN50'] = get_N50(lengths['gaps'], metrics['gapbases'])
            metrics['gapN90'] = get_N90(lengths['gaps'], metrics['gapbases'])
            metrics['meangap'] = get_mean(lengths['gaps'])
    if lengths['scaffolds']:
        metrics['maxscaffold'] = lengths['scaffolds'][-1]
        metrics['minscaffold'] = lengths['scaffolds'][0]
        if metrics['scaffoldbases']:  # scaffold N50
            metrics['scaffoldN50'] = get_N50(lengths['scaffolds'], 
                                             metrics['scaffoldbases'])
            metrics['scaffoldN90'] = get_N90(lengths['scaffolds'], 
                                           metrics['scaffoldbases'])
            metrics['meanscaffold'] = get_mean(lengths['scaffolds'])
    metrics['allbases'] = bases['total']
    metrics['maxlen'] = lengths['total'][-1]
    metrics['minlen'] = lengths['total'][0]
    metrics['N50'] = get_N50(lengths['total'], bases['total'])  # N50 of all


def basic_fasta_stats(fasta, min_gap, classic):
    '''Main method for stats calculation.  Creates data structures

       and controls workflow
    '''
    if not fasta:  # Assume STDIN
        fasta = sys.stdin
    else:
        fasta = return_filehandle(fasta)
    bases = {'A' : 0, 'a' : 0, 'C' : 0, 'c' : 0,
             'T' : 0, 't' : 0, 'G' : 0, 'g' : 0,
             'N' : 0, 'n' : 0, 'IUPAC' : 0, 'total' : 0}
    metrics = {'N50' : 0, 'maxlen': 0, 'minlen': 0,
               'contigN50' : 0, 'scaffoldN50' : 0,
               'contigN90': 0, 'scaffoldN90': 0, 'gapN90': 0,
               'meancontig': 0, 'meanscaffold': 0, 'meangap': 0,
               'gaps' : 0, 'gapN50' : 0, 'maxgap' : 0, 'mingap' : 0,
               'contigs' : 0, 'scaffolds' : 0, 'records' : 0,
               'maxscaffold' : 0, 'minscaffold' : 0,
               'maxcontig' : 0, 'mincontig' : 0,
               'contigbases' : 0, 'scaffoldbases' : 0, 'gapbases' : 0,
               'allbases' : 0, 'pgc' : 0}
    lengths = {'scaffolds' : [], 'gaps' : [], 'contigs' : [], 'total' : []}
    scheck = 0
    gcheck = 0
    ccheck = 0
    glen = 0
    clen = 0
    slen = 0
    length = 0
    i = 0
    for record in get_seqio_fasta_record(fasta):  # get records from SeqIO
        metrics['records'] += 1  # increment total
        seq = record.seq
        if seq:
            seq = seq.upper()
            length = len(seq)
            bases['total'] += length
            lengths['total'].append(length)
            while True:
                i = seq.find('N', i)  # split all sequences into contigs
                if i < 0:
                    if not scheck:  # not a scaffold yet
                        lengths['contigs'].append(length)
                        if length > 0:
#                            print(">{}\t{}\t{}".format(record.id, i, i))
                            metrics['contigs'] += 1
                            metrics['contigbases'] += length
                        for b in seq:
                            if b in bases:
                                bases[b] += 1
                            else:
                                bases['IUPAC'] += 1
                    else:  # scaffold
                        clen = length - ccheck
                        if clen > 0:
#                            print(">{}\t{}\t{}".format(record.id, i, i))
                            lengths['contigs'].append(clen)
                            metrics['contigs'] += 1
                            metrics['contigbases'] += clen
                        for b in seq[ccheck:]:
                            if b in bases:
                                bases[b] += 1
                            else:
                                bases['IUPAC'] += 1
                    break
                gcheck = i  # set position of gap start
                while i < length and seq[i] == 'N':  # get gap lengths
                    i += 1
                glen = i - gcheck
                if glen >= min_gap:  # check to see if this is a "gap"
                    scheck = 1  # set check for scaffold
                    clen = gcheck - ccheck
                    if clen > 0:
#                        print(">{}\t{}\t{}".format(record.id, i, i))
                        lengths['contigs'].append(clen)
                        metrics['contigs'] += 1
                        metrics['contigbases'] += clen
                    lengths['gaps'].append(glen)
                    bases['N'] += glen
                    metrics['gaps'] += 1
                    metrics['gapbases'] += glen
                    for b in seq[ccheck:gcheck]:
                        if b in bases:
                            bases[b] += 1
                        else:
                            bases['IUPAC'] += 1
                    ccheck = i
            if scheck:  # if scaffold as set by gap checks above
                metrics['scaffolds'] += 1
                lengths['scaffolds'].append(length)
                metrics['scaffoldbases'] += length
        scheck = 0  # reset all for the next sequence yielded
        gcheck = 0
        ccheck = 0
        glen = 0
        clen = 0
        i = 0
        pos = 0
    compile_metrics(metrics, lengths, bases)
    metrics['pgc'] = round((float(bases['G'] + bases['C'])/float(bases['total']))*100)
    if classic:
        metrics['scaffoldN50'] = metrics['N50']
        metrics['scaffolds'] = metrics['records']
        metrics['scaffoldbases'] = metrics['allbases']
        metrics['maxscaffold'] = metrics['maxlen']
        metrics['minscaffold'] = metrics['minlen']
        classic_metrics = OrderedDict([('Scaffolds', metrics['scaffolds']),
                            ('Max Scaffold', metrics['maxscaffold']),
                            ('Min Scaffold', metrics['minscaffold']),
                            ('Mean Scaffold', metrics['meanscaffold']),
                            ('Scaffold N50', metrics['scaffoldN50']),
                            ('Scaffold N90', metrics['scaffoldN90']),
                            ('Total Scaffold Length', metrics['scaffoldbases']),
                            ('Contigs', metrics['contigs']),
                            ('Max Contig', metrics['maxcontig']),
                            ('Min Contig', metrics['mincontig']),
                            ('Mean Contig', metrics['meancontig']),
                            ('Contig N50', metrics['contigN50']),
                            ('Contig N90', metrics['contigN90']),
                            ('Total Contig Length', metrics['contigbases']),
                            ('Assembly GC', metrics['pgc']),
                            ('Captured Gaps', metrics['gaps']),
                            ('Max Gap', metrics['maxgap']),
                            ('Min Gap', metrics['mingap']),
                            ('Mean Gap', metrics['meangap']),
                            ('Gap N50', metrics['gapN50']),
                            ('Total Gap Length', metrics['gapbases'])])
        return classic_metrics  # return GAEMR like output keys
    return metrics  # standard


@click.command()
@click.option('--fasta', help='''FASTA file to filter, can be compressed''')
@click.option('--classic', is_flag=True,
         help='''Outputs Stats with Older GAEMR Like Keys''')
@click.option('--human_readable', is_flag=True,
         help='''Outputs Human Readable Stats''')
@click.option('--min_gap', default=10, help="""Minimum length of consecutive N's to consider a gap and create a scaffold (default: 10)""")
@click.option('--log_file', default='./basic_fasta_stats.log',
help='''File to write log to.  (default:./basic_fasta_stats.log)''')
@click.option('--log_level', default='INFO',
help='''Log level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default:INFO)''')
def main(fasta, min_gap, classic, human_readable, log_file, log_level):
    '''Basic FASTA Stats Generation.  MORE DOC COMING'''
    log_level = getattr(logging, log_level.upper(), logging.INFO)
    msg_format = '%(asctime)s|%(name)s|[%(levelname)s]: %(message)s'
    logging.basicConfig(format=msg_format, datefmt='%m-%d %H:%M',
                        level=log_level)
    log_handler = logging.FileHandler(log_file, mode='w')
    formatter = logging.Formatter(msg_format)
    log_handler.setFormatter(formatter)
    logger = logging.getLogger('basic_fasta_stats')
    logger.addHandler(log_handler)
    if fasta and check_stdin(sys.stdin):
        logger.warning('stdin seen with FASTA, will process FASTA')
    if fasta:
        fasta = os.path.abspath(fasta)
    if human_readable:
        stats = basic_fasta_stats(fasta, min_gap, classic)
        for s in stats:
            print('{}\t{}'.format(s, stats[s]))
    else:
        print(json.dumps(basic_fasta_stats(fasta, min_gap, classic)))


if __name__ == '__main__':
    main()
