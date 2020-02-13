

import os
import sys
import re
import click
import json
import logging
from collections import OrderedDict
from time import sleep
from signal import signal, SIGPIPE, SIG_DFL
from ..helpers.sequence_helpers import get_seqio_fastq_record
from ..helpers.file_helpers import return_filehandle, check_stdin

signal(SIGPIPE, SIG_DFL) 


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


def compile_metrics(metrics, lengths, bases, bins, passes):
    '''Fill the metrics dictionary with the results from lengths and bases'''
    lengths['total'] = sorted(lengths['total'])
    metrics['allbases'] = bases['total']
    metrics['maxlen'] = lengths['total'][-1]
    metrics['minlen'] = lengths['total'][0]
    metrics['mean'] = get_mean(lengths['total'])  # N50 of all
    metrics['mean_passes'] = get_mean(lengths['passes'])
    metrics['passes_by_len'] = lengths['passes_by_len']  # redundant
    metrics['length_bins'] = [ (i, bins[i]) for i in sorted(bins.keys(),
                                                       key=lambda k: int(k)) ]
    metrics['passes_bins'] = [ (i, passes[i]) for i in sorted(passes.keys(), 
                                                        key=lambda k: int(k)) ]


def hifi_profiler(fastq, bin_size, split_passes):
    '''Main method for stats calculation.  Creates data structures

       and controls workflow
    '''
    if not fastq:  # Assume STDIN
        fastq = sys.stdin
    else:
        fastq = return_filehandle(fastq)
    bases = {'A': 0, 'a': 0, 'C': 0, 'c': 0,
             'T': 0, 't': 0, 'G': 0, 'g': 0,
             'N': 0, 'n': 0, 'IUPAC': 0, 'total': 0}
    metrics = {'maxlen': 0, 'minlen': 0, 'records': 0, 'mean_passes': 0,
               'mean': 0, 'allbases': 0}
    lengths = {'length_bins': [], 'passes': [], 'total': [], 
               'passes_bins': [], 'passes_by_len': {}}
    length = 0
    bins = {}  # for plotting
    passes = {}
    get_passes = re.compile('passes=(\d+)')
    for record in get_seqio_fastq_record(fastq):  # get records from SeqIO
        metrics['records'] += 1  # increment total
        seq = record.seq
        desc = record.description  # get header description to get passes
        if not desc:
            sys.stderr.write('NO DESCRIPTION FIELD NOT CALCULATING PASSES!\n')
        else:
            my_passes = get_passes.search(desc)
            my_passes = my_passes.groups(1)[0]
            lengths['passes'].append(int(my_passes))
            if my_passes not in passes:
                passes[my_passes] = 0
                lengths['passes_by_len'][my_passes] = []
            passes[my_passes] += 1
            lengths['passes_by_len'][my_passes].append(len(seq))
        if seq:
            seq = seq.upper()
            length = len(seq)  # cast as an int
            bases['total'] += length
            lengths['total'].append(length)
            bin_me = str(int(length/bin_size))  # bin number as an int
            if bin_me not in bins:
                bins[bin_me] = 0
            bins[bin_me] += 1
    compile_metrics(metrics, lengths, bases, bins, passes)
#    metrics['pgc'] = round((float(bases['G'] + bases['C'])/float(bases['total']))*100)
    return metrics  # standard


@click.command()
@click.option('--fastq', help='''FASTA file to filter, can be compressed''')
@click.option('--human_readable', is_flag=True,
         help='''Outputs Human Readable Stats''')
@click.option('--bin_size', default=1000, help="""Histogram Bin Size (default: 1000)""")
@click.option('--split_passes', is_flag=True, help="""Outputs reads into files based on the number of passes.""")
@click.option('--log_file', default='./hifi_profiler.log',
help='''File to write log to.  (default:./hifi_profiler.log)''')
@click.option('--log_level', default='INFO',
help='''Log level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default:INFO)''')
def main(fastq, human_readable, bin_size, split_passes, log_file, log_level):
    '''Reads HiFi data and produces metrics about passes.  MORE DOC COMING'''
    log_level = getattr(logging, log_level.upper(), logging.INFO)
    msg_format = '%(asctime)s|%(name)s|[%(levelname)s]: %(message)s'
    logging.basicConfig(format=msg_format, datefmt='%m-%d %H:%M',
                        level=log_level)
    log_handler = logging.FileHandler(log_file, mode='w')
    formatter = logging.Formatter(msg_format)
    log_handler.setFormatter(formatter)
    logger = logging.getLogger('hifi_profiler')
    logger.addHandler(log_handler)
    if fastq and check_stdin(sys.stdin):
        logger.warning('stdin seen with FASTQ, will process FASTQ')
    if fastq:
        fastq = os.path.abspath(fastq)
    stats = hifi_profiler(fastq, bin_size, split_passes)
    if human_readable:
        for s in stats:
            print('{}\t{}'.format(s, stats[s]))
    else:
        print(json.dumps(stats))


if __name__ == '__main__':
    main()
