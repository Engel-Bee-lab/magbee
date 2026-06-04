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
Check input contigs
"""
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

def get_sample_name_from_contig(filename):
        name = os.path.basename(filename)
        name = re.sub(r'\.(fa|fasta)(\.gz)?$', '', name)
        return name.split(".")[0]

if config['args']['mode'] == 'concatenate' and len(contig_files) > 1:
    raise ValueError("Multiple contig files found, but mode is set to 'concatenate'. Please provide a single concatenated contig file or change the mode to 'individual'.")
elif config['args']['mode'] == 'concatenate' and len(contig_files) == 1:
    concat_assembly = contig_files[0]

    # sample names will come from BAM files
    sample_names = []

elif config['args']['mode'] == 'individual' and len(contig_files) == 1:
    raise ValueError("Only one contig file found, but mode is set to 'individual'. Please provide separate contig files for each sample or change the mode to 'concatenate'.")

elif config['args']['mode'] == 'individual' and len(contig_files) > 1:
    contig_map = {}

    for f in contig_files:
        sample = get_sample_name_from_contig(f)
        if sample in contig_map:
            raise ValueError(f"Duplicate contig file for sample: {sample}")
        contig_map[sample] = f

    sample_names = list(contig_map.keys())
else:
    raise ValueError(f"Invalid mode: {config['args']['mode']}. Must be 'individual' or 'concatenate'.")

"""
Check bam files
"""
bam_dir = config['args']['bam_folder']

if not os.path.isdir(bam_dir):
    raise ValueError(f"BAM directory not found: {bam_dir}")

bam_files = glob.glob(os.path.join(bam_dir, "*_bam")) + \
            glob.glob(os.path.join(bam_dir, "*cluster_50")) 

if not bam_files:
    raise ValueError("No BAM files found")

# checking all the samples have bam files generated 
bam_map = {}

if config['args']['mode'] == 'concatenate':
    for bam_dir_path in bam_files:

        name = os.path.basename(bam_dir_path)

        # sample1_bam -> sample1
        sample = re.sub(r'_bam$', '', name)

        bam_inside = glob.glob(os.path.join(bam_dir_path, "*.bam"))

        if len(bam_inside) == 0:
            raise ValueError(
                f"No BAM file found in {bam_dir_path}"
            )

        if len(bam_inside) > 1:
            raise ValueError(
                f"Multiple BAM files found in {bam_dir_path}"
            )

        bam_map[sample] = bam_inside[0]

    sample_names = sorted(bam_map.keys())

elif config['args']['mode'] == 'individual':
    for f in bam_files:
        name = os.path.basename(f)

        matched = None
        for sample in sample_names:
            # safer prefix match
            if name.startswith(f"{sample}_"):
                matched = sample
                break

        if not matched:
            raise ValueError(f"Could not match BAM to any sample: {f}")

        if matched in bam_map:
            raise ValueError(f"Duplicate BAM for sample: {matched}")

        bam_map[matched] = f
    
    missing_bam = [
        s for s in sample_names
        if s not in bam_map
    ]

    if missing_bam:
        raise ValueError(
            f"Missing BAMs for samples: {missing_bam}"
        )

"""
Declaring other directories
"""
dir_assembly = contigs_dir
dir_backmapping = bam_dir
dir_reports = os.path.join(dir_out, 'REPORTS')
dir_binning = os.path.join(dir_out, 'PROCESSING' ,'5_binning')

"""Rules"""
include: os.path.join("rules", "5.metabat2.smk")
#include: os.path.join("rules", "5.concoct.smk")
include: os.path.join("rules", "5.vamb.smk")
#include: os.path.join("rules", "5.semibin.smk")
#include: os.path.join("rules", "5.comebin.smk")

include: os.path.join("rules", "6.bins_quality_checkm2.smk")

"""Mark target rules"""
target_rules = []
def targetRule(fn):
    assert fn.__name__.startswith('__')
    target_rules.append(fn.__name__[2:])
    return fn

"""
Defining the targets dictionary
"""
targets ={'binning':[], 'binning_qual':[]}

if config['args']['mode'] == 'individual':
    for sample in sample_names:
        #binning nightmare targets
        targets['binning'].append(os.path.join(dir_binning, "{sample}_metabat2_bins", "done.txt").format(sample=sample))
        targets['binning'].append(os.path.join(dir_binning, "all_metabat2_bins", "done.txt"))
        targets['binning'].append(os.path.join(dir_binning, "{sample}_vamb_bins", "done.txt").format(sample=sample))
        targets['binning'].append(os.path.join(dir_binning, "all_vamb_bins", "done.txt"))

        #these are erroing out in buidling training models, so not including tme for now. Testig GPU support for semibin2 training and comebin
        #targets['binning'].append(os.path.join(dir_binning, "{sample}_semibin_bins", "done.txt").format(sample=sample))
        #targets['binning'].append(os.path.join(dir_binning, "all_semibin_bins", "done.txt"))
        #targets['binning'].append(os.path.join(dir_binning, "{sample}_comebin_bins", "done.txt").format(sample=sample))

        targets['binning_qual'].append(os.path.join(dir_reports, "checkm2_all", "CheckM2_Metabat2_quality_report.tsv"))
        targets['binning_qual'].append(os.path.join(dir_reports, "checkm2_all", "CheckM2_VAMB_quality_report.tsv"))

        #targets['binning'].append(os.path.join(dir_binning, "{sample}_concoct", "bins", "done.txt").format(sample=sample))
        #targets['binning'].append(os.path.join(dir_binning, "all_conoct_bins", "renamed.txt"))
        #targets['binning_qual'].append(os.path.join(dir_reports, "checkm2_all", "CheckM2_CONCOCT_quality_report.tsv"))
            
        #targets['binning'].append(os.path.join(dir_binning, "gtdbtk_output", "done.txt"))
        #targets['binning'].append(os.path.join(dir_reports, "gtdbtk_bac120_summary.tsv")),

elif config['args']['mode'] == 'concatenate':
    for sample in sample_names:
        targets['binning'].append(os.path.join(dir_binning, "vamb_bins_concat", "done.txt").format(sample=sample))
        targets['binning'].append(os.path.join(dir_binning, "metabat2_bins_concat", "done.txt").format(sample=sample))


@targetRule
rule all:
    input: 
        targets['binning'],
        targets['binning_qual']