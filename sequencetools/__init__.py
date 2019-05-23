# -*- coding: utf-8 -*-
"""sequencetools -- Tools for FASTX manipulation and statistics"""

import os
import sys
import argparse
import logging
import click
from datetime import datetime
import locale
from signal import signal, SIGPIPE, SIG_DFL
from .version import version as VERSION

#
# Start coverage
#
#coverage.process_startup()
# set locale so grouping works
for localename in ['en_US', 'en_US.utf8', 'English_United_States']:
    try:
        locale.setlocale(locale.LC_ALL, localename)
        break
    except locale.Error:
        continue

AUTHOR = 'Connor Cameron'
EMAIL = 'ctc@ncgr.org'
COPYRIGHT = """Copyright (C) 2019, NCGR. All rights reserved.
"""
PROGRAM_NAME = 'sequencetools'
PROJECT_HOME = 'https://github.com/ncgr/seqtools2'
DEFAULT_FILE_LOGLEVEL = logging.DEBUG
DEFAULT_STDERR_LOGLEVEL = logging.INFO
STARTTIME = datetime.now()
#
# global logger object
#
logger = logging.getLogger(PROGRAM_NAME)
#
# private context function
#
_ctx = click.get_current_context


#
# Classes
#
class CleanInfoFormatter(logging.Formatter):
    def __init__(self, fmt='%(levelname)s: %(message)s'):
        logging.Formatter.__init__(self, fmt)

    def format(self, record):
        if record.levelno == logging.INFO:
            return record.getMessage()
        return logging.Formatter.format(self, record)


#
# Helper functions
#
#def init_dual_logger(file_log_level=DEFAULT_FILE_LOGLEVEL,
#                     stderr_log_level=DEFAULT_STDERR_LOGLEVEL):
#    '''Log to stderr and to a log file at different levels
#    '''
#    def decorator(f):
#        @functools.wraps(f)
#        def wrapper(*args, **kwargs):
#            global logger
#            # find out the verbose/quiet level
#            if _ctx().params['verbose']:
#                _log_level = logging.DEBUG
#            elif _ctx().params['quiet']:
#                _log_level = logging.ERROR
#            else:
#                _log_level = stderr_log_level
#            logger.setLevel(file_log_level)
#            stderrHandler = logging.StreamHandler(sys.stderr)
#            stderrFormatter = CleanInfoFormatter()
#            stderrHandler.setFormatter(stderrFormatter)
#            stderrHandler.setLevel(_log_level)
#            logger.addHandler(stderrHandler)
#            if _ctx().params['log']:  # start a log file in LOG_PATH
#                logfile_path = LOG_PATH / (PROGRAM_NAME + '.log')
#                if not LOG_PATH.is_dir():  # create LOG_PATH
#                    try:
#                        logfile_path.parent.mkdir(mode=0o755, parents=True)
#                    except OSError:
#                        logger.error('Unable to create log directory "%s"',
#                                     logfile_path.parent)
#                        raise OSError
#                else:
#                    if logfile_path.exists():
#                        try:
#                            logfile_path.unlink()
#                        except OSError:
#                            logger.error('Unable to remove log file "%s"',
#                                         logfile_path)
#                            raise OSError
#                logfileHandler = logging.FileHandler(str(logfile_path))
#                logfileFormatter = logging.Formatter(
#                    '%(levelname)s: %(message)s')
#                logfileHandler.setFormatter(logfileFormatter)
#                logfileHandler.setLevel(file_log_level)
#                logger.addHandler(logfileHandler)
#            logger.debug('Command line: "%s"', ' '.join(sys.argv))
#            logger.debug('%s version %s', PROGRAM_NAME, VERSION)
#            logger.debug('Run started at %s', str(STARTTIME)[:-7])
#            return f(*args, **kwargs)
#        return wrapper
#    return decorator


#def init_user_context_obj(initial_obj=None):
#    '''Put info from global options into user context dictionary
#    '''
#    def decorator(f):
#        @functools.wraps(f)
#        def wrapper(*args, **kwargs):
#            global config_obj
#            if initial_obj is None:
#                _ctx().obj = {}
#            else:
#                _ctx().obj = initial_obj
#            ctx_dict = _ctx().obj
#            if _ctx().params['verbose']:
#                ctx_dict['logLevel'] = 'verbose'
#            elif _ctx().params['quiet']:
#                ctx_dict['logLevel'] = 'quiet'
#            else:
#                ctx_dict['logLevel'] = 'default'
#            for key in ['progress', 'first_n']:
#                ctx_dict[key] = _ctx().params[key]
#            return f(*args, **kwargs)
#        return wrapper
#    return decorator

def cli():
    '''Sequencetools: FASTX Sequence Manipulation and Statistics

       USAGE: sequencetools <TOOL> [options]

       Current Tools:

         basic_fasta_stats
         format_fasta
         chunk_fasta
         chunk_fastq
         fastx_converter
         filter_fasta_by_length
         get_fasta_by_id
         get_fastq_by_id
         subset_fastq

       Please run `sequencetools <TOOL> --help` for individual usage
    '''
    if not len(sys.argv) > 1:
        print(cli.__doc__)
        sys.exit(1)
    tool = sys.argv[1].lower()  # get tool
    sys.argv = [tool] + sys.argv[2:]  # make usage for applicaiton correct
    if tool == '--help' or tool =='-h':
        print(cli.__doc__)
        sys.exit(1)
    if tool == 'basic_fasta_stats':
        from .tools import basic_fasta_stats
        basic_fasta_stats.main()
    if tool == 'format_fasta':
        from .tools import format_fasta
        format_fasta.main()
    if tool == 'chunk_fastq':
        from .tools import chunk_fastq
        chunk_fastq.main()
    if tool == 'chunk_fasta':
        from .tools import chunk_fasta
        chunk_fasta.main()
    if tool == 'fastx_converter':
        from .tools import fastx_converter
        fastx_converter.main()
    if tool == 'filter_fasta_by_length':
        from .tools import filter_fasta_by_length
        filter_fasta_by_length.main()
    if tool == 'get_fasta_by_id':
        from .tools import get_fasta_by_id
        get_fasta_by_id.main()
    if tool == 'get_fastq_by_id':
        from .tools import get_fastq_by_id
        get_fastq_by_id.main()
    if tool == 'subset_fastq':
        from .tools import subset_fastq
        subset_fastq.main()
