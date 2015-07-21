#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

execfile('HG1G2tools/version.py')

setup(
    name = 'HG1G2 tools',
    version = __version__,
    description = 'Implementation of the H,G1,G2 magnitude system in Python',
    long_description = open('README.txt').read(),
    author = 'Olli Wilkman',
    author_email = 'olli.wilkman@iki.fi',
    url = 'https://wiki.helsinki.fi/display/PSR/HG1G2+tools',
    packages = [
        'HG1G2tools'
    ],
    classifiers = [
        'TODO: Add trove classifiers (http://pypi.python.org/pypi?%3Aaction=list_classifiers)'
    ],
    license = "The MIT License",
    install_requires = ["numpy"]
)
