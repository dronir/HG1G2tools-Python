#!/usr/bin/env python

from distutils.core import setup

execfile('HG1G2tools/version.py')

setup(
    name = 'HG1G2 tools',
    version = __version__,
    description = 'TODO: Enter a description',
    long_description = open('README.rst').read() + '\n\n' + open('HISTORY.rst').read(),
    author = 'Olli Wilkman',
    author_email = 'olli.wilkman@iki.fi',
    url = 'TODO: Enter an URL',
    packages = [
        'HG1G2tools'
    ],
    classifiers = [
        'TODO: Add trove classifiers (http://pypi.python.org/pypi?%3Aaction=list_classifiers)'
    ]
)
