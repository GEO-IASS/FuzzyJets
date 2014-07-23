#!/usr/bin/python

import os, time, subprocess

def safe_mkdir(p):
    if not os.path.isdir(p):
        os.makedirs(p)

def submit_fuzzy(mu, size, lw, n_events, i):
    logdir = cwd + 'logs/' + time_postfix + "_bsub_" + str(mu)
    safe_mkdir(logdir)

    scratchdir = "/scratch/" + time_postfix + "_" + str(i)
    safe_mkdir(scratchdir)
    lw_flag = "1" if lw else "0"
    outfile_name = scratchdir + "/temp_" + str(size) + "s_" + lw_flag + "_" + str(i) + ".root"
    size_s = str(1.0*size / 10)
    submit = ['bsub', '-q', queue,
              '-R', 'select[(!preempt&&rhel60&&cvmfs&&inet)]',
              '-o', logdir + "/" + str(i) + ".log",
              subfile, workdir, scratchdir, outdir,
              "./Fuzzy", "--NEvents", str(n_events), "--OutFile", outfile_name,
              "--NPV", str(mu),  "--Size", str(size_s),
              "--LearnWeights", lw_flag, "--Batch", "1"]
    subprocess.call(submit)

user = "chstan"
workdir = "/u/at/chstan/nfs/summer_2014/ForConrad/"

cwd = os.getcwd() + "/"
subfile = cwd + "_batchSingleSub.sh"

time_postfix = time.strftime('%Y_%m_%d_%Hh%Mm')

events_per_job = 500
n_jobs = 20
queue = 'xlong'

outdir = cwd + 'files/' + time_postfix
safe_mkdir(outdir)

NPVs = [0, 10, 20, 30]
sizes = [7, 8, 9, 10]
learnWeights = [True, False]

for current_mu in NPVs:
    for current_size in sizes:
        for current_lw in learnWeights:
            for current_job in range(n_jobs):
                submit_fuzzy(current_mu, current_size, current_lw, events_per_job, current_job)
