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
Check input read files (same approach as Assembly_Snakefile)
"""
input_dir = config['args']['input']
extn = config['args']['extn']

sample = []

if config['args']['sequencing'] == 'paired':

    pattern_r1 = config['args']['pattern_r1']
    pattern_r2 = config['args']['pattern_r2']

    r1_files = glob.glob(os.path.join(input_dir, f"*{pattern_r1}*.{extn}"))
    r2_files = glob.glob(os.path.join(input_dir, f"*{pattern_r2}*.{extn}"))

    if not r1_files or not r2_files:
        raise ValueError("No R1 or R2 files found.")

    sample_inputs = {}

    def extract_sample_name(filename, ext, pattern_r1, pattern_r2):
        name = os.path.basename(filename)

        if name.endswith(f".{ext}"):
            name = name[:-(len(ext) + 1)]

        if name.endswith(pattern_r1):
            sample_name = name.rsplit(pattern_r1, 1)[0]
        elif name.endswith(pattern_r2):
            sample_name = name.rsplit(pattern_r2, 1)[0]
        else:
            raise ValueError(f"File does not end with R1/R2 pattern: {filename}")

        return sample_name

    for f in r1_files:
        sample_name = extract_sample_name(f, extn, pattern_r1, pattern_r2)
        sample_inputs.setdefault(sample_name, {})["r1"] = f

    for f in r2_files:
        sample_name = extract_sample_name(f, extn, pattern_r1, pattern_r2)
        sample_inputs.setdefault(sample_name, {})["r2"] = f

    paired_samples = {}

    for sample_name, reads in sample_inputs.items():
        if "r1" in reads and "r2" in reads:
            paired_samples[sample_name] = reads
        else:
            raise ValueError(f"Missing pair for sample {sample_name}")

    config["sample_names"] = paired_samples
    sample = sorted(paired_samples.keys())

"""
Check bins and collect MAG names
"""
bins_dir = config['args']['bins']

if not bins_dir or not os.path.isdir(bins_dir):
    raise ValueError(f"Bins directory not found: {bins_dir}")

bin_files = (
    glob.glob(os.path.join(bins_dir, "*.fa")) +
    glob.glob(os.path.join(bins_dir, "*.fasta")) +
    glob.glob(os.path.join(bins_dir, "*.fa.gz")) +
    glob.glob(os.path.join(bins_dir, "*.fasta.gz"))
)

if not bin_files:
    raise ValueError(f"No bin files found in: {bins_dir}")

def extract_mag_name(bin_file):
    name = os.path.basename(bin_file)
    name = re.sub(r'\.(fa|fasta)(\.gz)?$', '', name)
    return name

mags = sorted({extract_mag_name(f) for f in bin_files})


dir_reports = os.path.join(dir_out, 'REPORTS')
dir_species = os.path.join(dir_out, 'PROCESSING' ,'6_species_variation')

"""Rules"""
include: os.path.join("rules", "8.drep.smk")
include: os.path.join("rules", "9.instrain.smk")

"""
Defining the targets dictionary
"""
targets ={'derep':[], 'speciesVar':[]}

targets['derep'].append(os.path.join(dir_species, "drep_dastool", "done.txt"))


@targetRule
rule all:
    input: 
        targets['derep'],
        targets['speciesVar']

@targetRule
rule derep:
    input:
        targets['derep']

@targetRule
rule sepcies_variance:
    input:
        targets['speciesVar']
