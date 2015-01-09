#!/usr/bin/python

import os, time, itertools
from os import listdir
from os.path import isfile, join, getsize
from shutil import copyfile

repair_path = "/u/at/chstan/nfs/summer_2014/ForConrad/files/2500evts_qcd_ydm_offset_test/2015_01_06_14h43m15s"
#repair_path = "/example/"

damage_limit = 500
len_family = 5

files_and_sizes = [ (getsize(join(repair_path, f)), join(repair_path, f))
                    for f in listdir(repair_path)
                    if ("temp" in f) and isfile(join(repair_path,f)) ]

files = [f for (s,f) in files_and_sizes]

families = set([f[:-7] for (s, f) in files_and_sizes])
damaged_files = [ (s, f) for (s, f) in files_and_sizes if (s < damage_limit)]
damaged_families = set([f[:-7] for (s, f) in damaged_files])

def get_matching_strings(s, univ):
    return [t for t in univ if s in t]

print len(families)
for family in families:
    members = get_matching_strings(family, files)
    if len(members) != len_family:
        print "Incomplete family " + family
        continue

for family in damaged_families:
    members = get_matching_strings(family, files)
    damaged_members = [f for f in members if getsize(f) < damage_limit]
    undamaged_members = [f for f in members if getsize(f) >= damage_limit]
    if len(undamaged_members) == 0:
        print "Could not repair family " + family
        continue
    ps = [(d, u,) for d, u in
          itertools.izip(damaged_members, itertools.cycle(undamaged_members))]
    for d, u in ps:
        copyfile(u, d)
