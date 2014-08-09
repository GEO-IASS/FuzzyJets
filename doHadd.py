#!/usr/bin/python

import sys, pickle, subprocess, string

cleanscr = []

with open(sys.argv[1], 'rb') as f:
    cleanscr = pickle.load(f)

for command in cleanscr:
    command_arr = string.split(command, ' ')
    command_arr[0] = 'hadd -f'
    subprocess.call(' '.join(command_arr), shell = True)

# should also force remaking of files (put this in creation script)
# and should iterate through all the files in the directory clean_scripts
