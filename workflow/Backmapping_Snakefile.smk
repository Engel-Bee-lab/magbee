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
Check input reads
"""
input_dir = config['args']['input']

# List of file paths matching the pattern
#replace R1 to 1 for SRA reads
extn=config['args']['extn']

if config['args']['sequencing'] == 'paired':

    pattern_r1 = config['args']['pattern_r1']
    pattern_r2 = config['args']['pattern_r2']

    # Step 1: Find files
    r1_files = glob.glob(os.path.join(input_dir, f"*{pattern_r1}.*{extn}"))
    r2_files = glob.glob(os.path.join(input_dir, f"*{pattern_r2}.*{extn}"))

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
Check input contigs
"""
##Here checking the contigs input, each sample should have a contigs folder here
contigs_dir = config['args']['contigs']

if not os.path.isdir(contigs_dir):
    raise ValueError(f"Contigs directory not found: {contigs_dir}")

# find contig files (.fa / .fasta)
contig_files = glob.glob(os.path.join(contigs_dir, "*.fa")) + \
               glob.glob(os.path.join(contigs_dir, "*.fasta")) + \
               glob.glob(os.path.join(contigs_dir, "*.fa.gz")) + \
               glob.glob(os.path.join(contigs_dir, "*.fasta.gz"))

if not contig_files:
    raise ValueError("No contig files found")

# map sample -> contig file
def get_sample_name_from_contig(filename):
    name = os.path.basename(filename)

    # remove extensions (.fa/.fasta + optional .gz)
    name = re.sub(r'\.(fa|fasta)(\.gz)?$', '', name)

    # take only the first token before a dot
    sample = name.split(".")[0]

    return sample

contig_map = {}

for f in contig_files:
    name = os.path.basename(f)
    name = re.sub(r'\.(fa|fasta)(\.gz)?$', '', name)

    if name == "concatenated_assemblies":
        continue

    matched = None
    for sample in sample_names:
        if name.startswith(sample):
            matched = sample
            break

    if not matched:
        raise ValueError(f"Could not match contig file to any sample: {f}")

    if matched in contig_map:
        raise ValueError(f"Duplicate contig file for sample: {matched}")

    contig_map[matched] = f

# validate: every sample has contigs
missing = [s for s in sample_names if s not in contig_map]

if missing:
    raise ValueError(f"Missing contigs for samples: {missing}")

# optional: store in config for later rules
config["contigs"] = contig_map

"""
Declaring directories for each step
"""
#making directories for each step
#Saving most of the files to PROCESSING, sine they are intermediate files
#dir_assembly = os.path.join(dir_out, 'PROCESSING' ,'3_coassembly')
dir_hostcleaned= input_dir
dir_assembly= contigs_dir
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

#############################################
# AFTER sample detection logic for maping 
#############################################

#sample_names = list(paired_samples.keys())
N_SAMPLES = len(sample_names)
THRESHOLD = config['args'].get('mapping_threshold', 50)

#print(f"Detected {N_SAMPLES} samples")

"""
Rules
"""
if config['args']['mode'] == 'individual':
    print("Using individual mode")
    if N_SAMPLES <= THRESHOLD:
        print("Using all-vs-all mapping strategy")
        include: os.path.join("rules", "4.mapping_all_2_all.smk")
    else:
        print("Using simka mapping strategy")
        include: os.path.join("rules", "4.mapping_simka.smk")
        include: os.path.join("rules", "4.backmapping_simka.smk")
elif config['args']['mode'] == 'concatenate':
    print("Using concatenate mode")
    include: os.path.join("rules", "4.backmapping_concat.smk")
include: os.path.join("rules", "Final_backmapping.smk")

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
    if config['args']['mode'] == 'individual':
        for sample in sample_names:
            if N_SAMPLES <= THRESHOLD:
                targets['backmapping'].append(os.path.join(dir_backmapping, "{sample}_index.mmi").format(sample=sample))
                targets['backmapping'].append(os.path.join(dir_backmapping, "{sample}_bam", "done.txt").format(sample=sample))
                targets['backmapping'].append(os.path.join(dir_reports, "backmapping_report_all_to_all.txt"))
            else:
                targets['backmapping'].append(os.path.join(dir_backmapping, "simka_input_list.txt"))
                targets['backmapping'].append(os.path.join(dir_backmapping, "{sample}_cluster_50", "{sample}_backmap_samples.txt").format(sample=sample))
                targets['backmapping'].append(os.path.join(dir_backmapping, "{sample}_cluster_50", "done.txt").format(sample=sample))
                targets['backmapping'].append(os.path.join(dir_reports, "backmapping_report_simka.txt"))
    elif config['args']['mode'] == 'concatenate':
        targets['backmapping'].append(os.path.join(dir_assembly, "concatenated_assemblies.fa.gz"))
        targets['backmapping'].append(os.path.join(dir_assembly, "concatenated_assemblies.mmi"))
        for sample in sample_names:
            targets['backmapping'].append(os.path.join(dir_backmapping, "{sample}_bam", "done.txt").format(sample=sample))
        targets['backmapping'].append(os.path.join(dir_reports, "backmapping_report_all_to_concat.txt"))
    else:
        raise ValueError(f"Invalid mode: {config['args']['mode']}")
else:
    raise ValueError(f"Invalid sequencing type: {config['args']['sequencing']}")

@targetRule
rule all:
    input:
        targets['backmapping']