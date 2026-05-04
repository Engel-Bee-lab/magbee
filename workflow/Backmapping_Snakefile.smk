import yaml
import os
import glob
import re
import sys
import shutil
from metasnek import fastq_finder

"""Parse config"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")
"""
PREFLIGHT CHECKS
"""
"""
Declaring directories
"""
dir = {}
#declaring output file
try:
    if config['args']['output'] is None:
        dir_out = os.path.join('output')
    else:
	    dir_out = config['args']['output']
except KeyError:
    dir_out = os.path.join('output')

# temp dir
if config['args']['temp_dir'] is None:
    dir_temp = os.path.join(dir_out, "temp")
else:
    dir_temp = config['args']['temp_dir']

#declaring some the base directories
dir_env = os.path.join(workflow.basedir,"envs")
dir_script = os.path.join(workflow.basedir,"scripts")

"""
Check input files
"""
input_dir = config['args']['contigs_dir']

fa_files = glob.glob(os.path.join(input_dir, "*.fa"))
fasta_files = glob.glob(os.path.join(input_dir, "*.fasta"))

seq_files = fa_files + fasta_files

if not seq_files:
    raise ValueError("No .fa or .fasta files found")

sample_inputs = {}

for f in seq_files:
    name = os.path.basename(f)
    sample = name.rsplit(".", 1)[0]

    if sample in sample_inputs:
        raise ValueError(f"Duplicate sample detected: {sample}")

    sample_inputs[sample] = f

config["sample_names"] = sample_inputs

"""
Declaring directories for each step
"""
#making directories for each step
#Saving most of the files to PROCESSING, sine they are intermediate files
#dir_assembly = os.path.join(dir_out, 'PROCESSING' ,'3_coassembly')
dir_assembly= input_dir
dir_backmapping = os.path.join(dir_out, 'PROCESSING' ,'4_backmapping')
dir_reports = os.path.join(dir_out, 'REPORTS')

"""ONSTART/END/ERROR
Tasks to perform at various stages the start and end of a run.
"""
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    target_log = os.path.join(dir['log'], os.path.basename(current_log))
    shutil.copy(current_log, target_log)

dir = {'log': os.path.join(dir_out, 'logs')}

def cleanup_logs():
    if os.path.isdir(dir["log"]):
        oldLogs = filter(re.compile(r'.*.log').match, os.listdir(dir["log"]))
        for logfile in oldLogs:
            os.unlink(os.path.join(dir["log"], logfile))

"""Rules"""
if N_SAMPLES <= THRESHOLD:
    print("Using all-vs-all mapping strategy"),
    include: os.path.join("rules", "4.mapping_all_2_all.smk")
else:
    print("Using simka mapping strategy"),
    include: os.path.join("rules", "4.mapping_simka.smk")
    include: os.path.join("rules", "4.backmapping_simka.smk")

#############################################
# AFTER sample detection logic for maping 
#############################################

#sample_names = list(paired_samples.keys())
N_SAMPLES = len(sample_names)
THRESHOLD = config['args'].get('mapping_threshold', 50)

#print(f"Detected {N_SAMPLES} samples")


"""Mark target rules"""
target_rules = []
def targetRule(fn):
    assert fn.__name__.startswith('__')
    target_rules.append(fn.__name__[2:])
    return fn

"""
Defining the targets dictionary
"""
targets ={'backmapping':[]}

if config['args']['sequencing'] == 'paired':
    for sample in sample_names:
        if N_SAMPLES <= THRESHOLD:
            targets['backmapping'].append(os.path.join(dir_backmapping, "{sample}_index.mmi").format(sample=sample))
            targets['backmapping'].append(os.path.join(dir_backmapping, "{sample}_bam", "done.txt").format(sample=sample))
        else:
            targets['backmapping'].append(os.path.join(dir_backmapping, "simka_input_list.txt"))
            targets['backmapping'].append(os.path.join(dir_backmapping, "{sample}_cluster_50", "{sample}_backmap_samples.txt").format(sample=sample))
            targets['backmapping'].append(os.path.join(dir_backmapping, "{sample}_cluster_50", "done.txt").format(sample=sample))

@targetRule
rule all:
    input:
        targets['backmapping']