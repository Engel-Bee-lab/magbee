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
input_dir = config['args']['input']

# List of file paths matching the pattern
#replace R1 to 1 for SRA reads
extn=config['args']['extn']

if config['args']['sequencing'] == 'paired':

    pattern_r1 = config['args']['pattern_r1']
    pattern_r2 = config['args']['pattern_r2']

    # Step 1: Find files
    r1_files = glob.glob(os.path.join(input_dir, f"*{pattern_r1}*.{extn}"))
    r2_files = glob.glob(os.path.join(input_dir, f"*{pattern_r2}*.{extn}"))

    if not r1_files or not r2_files:
        raise ValueError("No R1 or R2 files found.")

    # Step 2: Build mapping
    sample_inputs = {}

    def extract_sample_name(filename, ext, pattern_r1, pattern_r2):
        name = os.path.basename(filename)
        
        # remove extension
        if name.endswith(f".{ext}"):
            name = name[:-(len(ext) + 1)]
        
        # determine which pattern is at the end
        if name.endswith(pattern_r1):
            sample = name.rsplit(pattern_r1, 1)[0]
        elif name.endswith(pattern_r2):
            sample = name.rsplit(pattern_r2, 1)[0]
        else:
            raise ValueError(f"File does not end with R1/R2 pattern: {filename}")
        
        return sample
            
    for f in r1_files:
        sample = extract_sample_name(f, extn, pattern_r1, pattern_r2)
        sample_inputs.setdefault(sample, {})["r1"] = f

    for f in r2_files:
        sample = extract_sample_name(f, extn, pattern_r1, pattern_r2)
        sample_inputs.setdefault(sample, {})["r2"] = f

    # Step 3: Validate pairs
    paired_samples = {}

    for sample, reads in sample_inputs.items():
        print (sample, reads)
        if "r1" in reads and "r2" in reads:
            paired_samples[sample] = reads
        else:
            raise ValueError(f"Missing pair for sample {sample}")

    config["sample_names"] = paired_samples
    sample_names = list(paired_samples.keys())
    #print(f"Sample inputs: {paired_samples}")

"""
Declaring directories for each step
"""
#making directories for each step
#Saving most of the files to PROCESSING, sine they are intermediate files
dir_fastp = os.path.join(dir_out, 'PROCESSING' ,'1_fastp')
dir_hostcleaned = os.path.join(dir_out, 'PROCESSING' ,'2_host_cleaned')
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

onstart:
    cleanup_logs()

onsuccess:
    sys.stderr.write('\n\nWorkflow ran successfully!\n\n')

onerror:
    sys.stderr.write('\n\nWorkflow run failed\n\n')

"""Rules"""
include: os.path.join("rules", "1.qc.smk")
include: os.path.join("rules", "2.host_contamination.smk")
include: os.path.join("rules", "Final_QC.smk")

"""Mark target rules"""
target_rules = []
def targetRule(fn):
    assert fn.__name__.startswith('__')
    target_rules.append(fn.__name__[2:])
    return fn

"""
Defining the targets dictionary
"""
targets ={'qc':[]}

if config['args']['sequencing'] == 'paired':
    for sample in sample_names:
        targets['qc'].append(expand(os.path.join(dir_fastp,"{sample}.stats.html"), sample=sample)),
        targets['qc'].append(expand(os.path.join(dir_hostcleaned,"{sample}_temp.bam"), sample=sample)),
        targets['qc'].append(expand(os.path.join(dir_hostcleaned,"{sample}_bamstats.txt"), sample=sample)),
        targets['qc'].append(expand(os.path.join(dir_hostcleaned,"{sample}_R1.hostcleaned.fastq.gz"), sample=sample)),
        targets['qc'].append(expand(os.path.join(dir_hostcleaned,"{sample}_R2.hostcleaned.fastq.gz"), sample=sample)),
        targets['qc'].append(expand(os.path.join(dir_reports,"QC","{sample}_host_qc_report.txt"), sample=sample)),
        targets['qc'].append(os.path.join(dir_reports, "final_host_qc_summary.txt"))

@targetRule
rule all:
    input:
        targets['qc']

