#!/bin/python3
# Python master file to parse args and call the de-multiplexer
import argparse
import subprocess

__author__ = "Jonathan Smith"

# 1) Handle input arguments
parser = argparse.ArgumentParser(description="Process input arguments for Cassava.py.")
parser.add_argument("-test", "-test_flag", "--test", "--test_flag")
args = parser.parse_args()

# 2) Validate input arguments - only in exceptional flags

# 3) Submit cassava job
command = \
    '''
        ssh 10.202.64.28 "
        source /etc/profile;
        qsub
            -v testflag='%s'
            test-script.sh
    ''' \
        % (args.get("test"))
command = ' '.join(command.split())

process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE).stdout
out = process.readlines()

# 4) Setup monitoring and waiting

# 5) File our outputs

# 6) File our displays

# 7) Cleanup the environment
