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
