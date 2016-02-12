#!/usr/bin/env bash

# Executable sh for cassava based de-multiplexing
# Assumptions that are made include:
#   all necessary files are on this path
#   all necessary arguments are present

# setup python environment
module add python3

# Pass direct to python as its easier to parse and validate there - no need to qsub though
python cassava.py "$@"