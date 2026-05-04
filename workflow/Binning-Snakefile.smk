import yaml
import os
import glob

"""Parse config"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")

"""Rules"""
include: os.path.join("rules", "5.metabat2.smk")
include: os.path.join("rules", "5.concoct.smk")
include: os.path.join("rules", "5.vamb.smk")
#include: os.path.join("rules", "5.semibin.smk")
#include: os.path.join("rules", "5.comebin.smk")

include: os.path.join("rules", "6.bins_quality.smk")

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


if config['args']['sequencing'] == 'paired':
    for sample in sample_names:
        #binning nightmare targets
        targets['binning'].append(os.path.join(dir_binning, "{sample}_metabat2_bins", "done.txt").format(sample=sample))
        targets['binning'].append(os.path.join(dir_binning, "all_metabat2_bins", "done.txt"))
        targets['binning'].append(os.path.join(dir_binning, "{sample}_vamb_bins", "done.txt").format(sample=sample))
        targets['binning'].append(os.path.join(dir_binning, "all_vamb_bins", "done.txt"))
        targets['binning'].append(os.path.join(dir_binning, "{sample}_concoct", "bins", "done.txt").format(sample=sample))
        targets['binning'].append(os.path.join(dir_binning, "all_conoct_bins", "renamed.txt"))
        
        #these are erroing out in buidling training models, so not including tme for now.
        #targets['binning'].append(os.path.join(dir_binning, "{sample}_semibin_bins", "done.txt").format(sample=sample))
        #targets['binning'].append(os.path.join(dir_binning, "all_semibin_bins", "done.txt"))
        #targets['binning'].append(os.path.join(dir_binning, "{sample}_comebin_bins", "done.txt").format(sample=sample))

        targets['binning_qual'].append(os.path.join(dir_binning, "checkm2", "checkm2_output_metabat2", "quality_report.tsv"))
        targets['binning_qual'].append(os.path.join(dir_binning, "checkm2", "checkm2_output_concoct", "quality_report.tsv"))
        targets['binning_qual'].append(os.path.join(dir_binning, "checkm2", "checkm2_output_vamb", "quality_report.tsv"))
        
        #targets['binning'].append(os.path.join(dir_binning, "gtdbtk_output", "done.txt"))
        #targets['binning'].append(os.path.join(dir_reports, "gtdbtk_bac120_summary.tsv")),

@targetRule
rule all:
    input: 
        targets['binning'],
        targets['binning_qual']