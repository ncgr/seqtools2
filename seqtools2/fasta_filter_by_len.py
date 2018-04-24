#!/usr/bin/env python

import sys
from Bio import SeqIO

parser = argparse.ArgumentParser(description='''
''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--genome', metavar = '</path/to/my/genome.fna>',
help='''Genome to check''')

parser.add_argument('--annotation', metavar = '</path/to/my/annotation.gff>',
help='''Annotation to check''')


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
logger = logging.getLogger('detect_incongruencies')
logger.addHandler(log_handler)


if __name__ == '__main__':
    
