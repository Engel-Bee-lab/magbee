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
dir_hostcleaned = config['args']['input']
dir_reports = os.path.join(dir_out, 'REPORTS')
dir_binning = os.path.join(dir_out, 'PROCESSING' ,'5_binning')
dir_species = os.path.join(dir_out, 'PROCESSING' ,'6_species_variation')

"""
Check input read files (same approach as Assembly_Snakefile)
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

        # remove common read suffixes before matching the read pattern
        suffixes = [f".{ext}", ".fastq", ".fq", ".hostcleaned", ".trimmed"]
        changed = True
        while changed:
            changed = False
            for suffix in suffixes:
                if name.endswith(suffix):
                    name = name[: -len(suffix)]
                    changed = True
                    break
        
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
        #print ("samples and reads found")
        #print (sample, reads)
        if "r1" in reads and "r2" in reads:
            paired_samples[sample] = reads
        else:
            raise ValueError(f"Missing pair for sample {sample}")

    config["sample_names"] = paired_samples
    sample_names = list(paired_samples.keys())
    #print(f"[DEBUG] sample names: {sample_names}")
    #print(f"Sample inputs: {paired_samples}")

"""
Check bins and collect MAG names
"""
bins_dir = config['args']['bins']


"""Rules"""
include: os.path.join("rules", "8.drep.smk")
include: os.path.join("rules", "9.instrain.smk")
#include: os.path.join("rules", "10.bakta.smk")

"""Mark target rules"""
target_rules = []

def targetRule(fn):
    assert fn.__name__.startswith('__')
    target_rules.append(fn.__name__[2:])
    return fn

"""
Defining the targets dictionary
"""
targets ={'derep':[], 'speciesVar':[]}
targets['derep'].append(os.path.join(dir_species, "drep_dastools", "done.txt"))
#targets['speciesVar'].extend(
#    expand(
#        os.path.join(dir_species, "bakta", "{derep_genome}.done"),
#        derep_genome=get_derep_bin_names(),
#    )
#)
targets['speciesVar'].extend(
    expand(
        os.path.join(dir_species, "inStrain", "instrain_profile_db_mode", "{sample}_profile_db_mode.done"),
        sample=sample_names,
    )
)
if sample_names:
    targets['speciesVar'].append(os.path.join(dir_species, "inStrain", "instrain_compare", "all_compared.done"))

targets['speciesVar'].append(os.path.join(dir_species, "inStrain", "relative_abundance", "All_genome_info.tsv")),
targets['speciesVar'].append(os.path.join(dir_species, "inStrain", "relative_abundance", "rel_abundance_long.tsv")),
targets['speciesVar'].append(os.path.join(dir_species, "inStrain", "relative_abundance", "rel_abundance_matrix.tsv")),
targets['speciesVar'].append(os.path.join(dir_species, "inStrain", "relative_abundance", "prevalence.tsv"))

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
