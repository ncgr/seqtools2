# -*- coding: utf-8 -*-
#
"""sequencetools -- A package for FASTX stats generation and manipulation"""
#
# Developers, install with:
#    pip install -r requirements.txt
#    python setup.py develop
#
import sys
from setuptools import setup, find_packages
# Version restrictions and dependencies
if sys.version_info < (3, 4, 0, 'final', 0):
    raise SystemExit("This package requires python 3.4 or higher.")

NAME = 'sequencetools'

tests_require = []
#tests_require = [
#    'check-manifest>=0.25',
#    'isort>=4.2.2.2',
#    'pydocstyle>=1.0.0',
#    'pytest-cache>=1.0',
#    'pytest-cov>=1.8.0',
#    'pytest-pep8>=1.0.6',
#    'pytest>=2.8.0'
#]

extras_require = dict(docs=['Sphinx>=1.4.2'], tests=tests_require)

extras_require['all'] = []
#for reqs in extras_require.values():
#    extras_require['all'].extend(reqs)

packages = find_packages()

setup(
    description=__doc__,
    packages=packages,
    setup_requires=['setuptools>30.3.0',
                    'setuptools-scm>1.5'
                    ],
    entry_points={
        'console_scripts': [NAME + ' = ' + NAME + ':cli']
    },
    version = '0.0.3',
#    use_scm_version={
#        'version_scheme': 'guess-next-dev',
#        'local_scheme': 'dirty-tag',
#        'write_to': NAME + '/version.py'
#    },
    extras_require=extras_require,
    tests_require=tests_require,
)
