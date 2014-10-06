#!/usr/bin/python

import os, time, subprocess, pickle, threading

def safe_mkdir(p):
    if not os.path.isdir(p):
        os.makedirs(p)

def file_name(prefix, is_intermediate, size, lw, mu, rec, cut):
    lw_flag = "1" if lw else "0"
    rec_flag = "1" if rec else "0"
    div = "/temp_" if is_intermediate else "/"
    return prefix + div + str(size) + "s_" + str(mu) + "mu_" + lw_flag + "lw_" + rec_flag + "rec_" + str(cut) + "cut"

def subprocess_call_bg(args):
    subprocess.call(args)

def submit_fuzzy(mu, size, lw, rec, cut, n_events, i, unique_id):
    logdir = cwd + 'logs/' + name + "/" + time_postfix + "_bsub_" + str(mu) + "_" + str(unique_id)
    safe_mkdir(logdir)

    scratchdir = "/scratch/chstan/" + time_postfix + "_" + str(unique_id)
    safe_mkdir(scratchdir)
    lw_flag = "1" if lw else "0"
    rec_flag = "1" if rec else "0"
    outfile_name = file_name(scratchdir, True, size, lw, mu, rec, cut) + "_" + str(i) + ".root"
    size_s = str(1.0*size / 10)
    submit = ['bsub', '-q', queue,
              '-R', 'select[(!preempt&&rhel60&&cvmfs&&inet)]',
              '-o', logdir + "/" + str(i) + ".log",
              subfile, workdir, scratchdir, outdir,
              "./Fuzzy", "--NEvents", str(n_events), "--OutFile", outfile_name,
              "--NPV", str(mu),  "--Size", str(size_s),
              "--PythiaConfig", pythia_conf,
              "--LearnWeights", lw_flag, "--Recombination", rec_flag,
              "--pTMin", str(cut), "--Batch", "1"]
    subprocess.call(submit)

workdir = "/u/at/chstan/nfs/summer_2014/ForConrad/"

pythia_conf = workdir + "configs/default_batch.pythia"

cwd = os.getcwd() + "/"
subfile = cwd + "_batchSingleSub.sh"

time_postfix = time.strftime('%Y_%m_%d_%Hh%Mm%Ss')

events_per_job = 5000
n_jobs = 50
n_jobs_patch = 10
queue = 'xlong'

name = '250kevts_zprime_mu0'
outdir = cwd + 'files/' + name + '/' + time_postfix
safe_mkdir(outdir)

NPVs = [0]
sizes = [10]
learnWeights = [False, True]
recombinations = [False, True]
pT_cuts = [5, 15, 25, 50]

j = 0
cleanup_commands = []

for current_mu in NPVs:
    for current_size in sizes:
        for current_lw in learnWeights:
            for current_rec in recombinations:
                for current_pT_cut in pT_cuts:
                    out_base = file_name(outdir, False, current_size, current_lw, current_mu, current_rec, current_pT_cut)
                    cleanup_base = file_name(outdir, True, current_size, current_lw, current_mu, current_rec, current_pT_cut)
                    cleanup_command = "hadd " + out_base + ".root" + " " \
                        + cleanup_base + "_{0.." + str(n_jobs-1) + "}.root"
                    cleanup_commands.append(cleanup_command)
                    for current_job in range(n_jobs):
                        submit_fuzzy(current_mu, current_size, current_lw, current_rec, current_pT_cut, events_per_job, current_job, j)
                        j += 1

outdir = cwd + 'files/' + name + '_patch/' + time_postfix
safe_mkdir(outdir)
for current_mu in NPVs:
    for current_size in sizes:
        for current_lw in learnWeights:
            for current_rec in recombinations:
                for current_pT_cut in pT_cuts:
                    for current_job in range(n_jobs_patch):
                        submit_fuzzy(current_mu, current_size, current_lw, current_rec, current_pT_cut, events_per_job, current_job, j)
                        j += 1

with open('clean_scripts/' + time_postfix + '.clscr', 'wb') as outf:
    pickle.dump(cleanup_commands, outf)

print "SUBMITTED " + str(len(NPVs) * len(learnWeights) * len(sizes) * (n_jobs + n_jobs_patch) * len(recombinations) * len(pT_cuts)) + " JOBS TO THE QUEUE " + queue
