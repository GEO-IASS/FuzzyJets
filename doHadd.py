#!/usr/bin/python

import sys, pickle, subprocess

with open(sys.argv[1], 'rb') as f:
    cleanscr = pickle.load(f)

for command in cleanscr:
    subprocess.call(command, shell = True)
