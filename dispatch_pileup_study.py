#!/usr/bin/python

import os, time, subprocess, pickle, threading, random

def safe_mkdir(p):
    if not os.path.isdir(p):
        os.makedirs(p)

def file_name(prefix, is_intermediate):
    div = "/temp" if is_intermediate else "/fin"
    return prefix + div

def submit_fuzzy(i, unique_id):
    logdir = cwd + 'logs/' + name + "/" + time_postfix + "_bsub_" + str(unique_id)
    safe_mkdir(logdir)

    scratchdir = "/scratch/chstan/" + time_postfix + "_" + str(unique_id)
    safe_mkdir(scratchdir)
    outfile_name = file_name(scratchdir, True) + "_" + str(i) + ".root"
    submit = ['bsub', '-q', queue,
              '-R', 'select[(!preempt&&rhel60&&cvmfs&&inet)]',
              '-o', logdir + "/" + str(i) + ".log",
              subfile, workdir, scratchdir, outdir,
              "./Fuzzy", "--OutFile", outfile_name,
              "--Size", str(1),
              "--NPV", str(0),
              "--PythiaConfig", pythia_conf,
              "--Batch", "1",
              "--SigmaDependenceOnMuCompleteStudy", "1"]
    #print " ".join(submit)
    time.sleep(random.randint(delay_secs_low, delay_secs_high))
    subprocess.call(submit)

workdir = "/u/at/chstan/nfs/summer_2014/ForConrad/"

pythia_conf = workdir + "configs/default_batch.pythia"

cwd = os.getcwd() + "/"
subfile = cwd + "_batchSingleSub.sh"

time_postfix = time.strftime('%Y_%m_%d_%Hh%Mm%Ss')

n_jobs = 100
n_jobs_patch = 5
delay_secs_low = 0
delay_secs_high = 0
queue = 'xlong'

name = '20k_zprime_pileup_effect_corrected'
outdir = cwd + 'files/' + name + '/' + time_postfix
safe_mkdir(outdir)


j = 0
cleanup_commands = []


out_base = file_name(outdir, False)
cleanup_base = file_name(outdir, True)
cleanup_command = "hadd " + out_base + ".root" + " " \
                  + cleanup_base + "_{0.." + str(n_jobs-1) + "}.root"
cleanup_commands.append(cleanup_command)
for current_job in range(n_jobs):
    submit_fuzzy(current_job, j)
    j += 1

outdir = cwd + 'files/' + name + '_patch/' + time_postfix
safe_mkdir(outdir)

for current_job in range(n_jobs_patch):
    submit_fuzzy(current_job, j)
    j += 1

with open('clean_scripts/' + time_postfix + '.clscr', 'wb') as outf:
    pickle.dump(cleanup_commands, outf)

print "SUBMITTED " + str((n_jobs + n_jobs_patch)) + " JOBS TO THE QUEUE " + queue
