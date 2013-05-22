#!/usr/bin/env bash
echo $1 > VERSION
echo "__version__ = '$1'" > HG1G2tools/version.py
