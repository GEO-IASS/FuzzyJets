#!/usr/bin/python

import sys, pickle, subprocess

with open(sys.argv[1], 'rb') as f:
    cleanscr = pickle.load(f)

for command in cleanscr:
    subprocess.call(command, shell = True)

# should also force remaking of files (put this in creation script)
# and should iterate through all the files in the directory clean_scripts
