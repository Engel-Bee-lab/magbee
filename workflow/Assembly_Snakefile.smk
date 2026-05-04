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
#dir_assembly = os.path.join(dir_out, 'PROCESSING' ,'3_coassembly')
dir_hostcleaned= input_dir
dir_assembly = os.path.join(dir_out, 'PROCESSING' ,'3_individual_assembly')
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
include: os.path.join("rules", "3.individual_assembly.smk")
include: os.path.join("rules", "Final_assembly.smk")

"""Mark target rules"""
target_rules = []
def targetRule(fn):
    assert fn.__name__.startswith('__')
    target_rules.append(fn.__name__[2:])
    return fn

"""
Defining the targets dictionary
"""
targets ={'assemble':[]}

if config['args']['sequencing'] == 'paired':
    ## For individual assembly targets
    for sample in sample_names:
        targets['assemble'].append(os.path.join(dir_assembly,"{sample}.megahit.contigs.fa").format(sample=sample)),
        targets['assemble'].append(os.path.join(dir_assembly, "{sample}_quast_output", "report.txt").format(sample=sample)),
        targets['assemble'].append(os.path.join(dir_assembly, "{sample}_assembly_report.txt").format(sample=sample)),
        targets['assemble'].append(os.path.join(dir_reports, "Assembly_stats_all.csv")),
        targets['assemble'].append(os.path.join(dir_reports, "assembly", "{sample}.megahit.contigs.fa.gz").format(sample=sample))

@targetRule
rule all:
    input:
        targets['assemble']