import yaml
import os
import glob

"""Parse config"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")

"""Rules"""
#include: os.path.join("rules", "3.co-assembly.smk")
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
    ## For coassembly targets
    #targets['assemble'].append(os.path.join(dir_assembly,"coassembly.megahit.contigs.fa"))
    #targets['assemble'].append(os.path.join(dir_assembly,"quast_report.txt"))

    ## For individual assembly targets
    for sample in sample_names:
        targets['assemble'].append(os.path.join(dir_assembly,"{sample}.megahit.contigs.fa").format(sample=sample)),
        targets['assemble'].append(os.path.join(dir_assembly, "{sample}_quast_output", "report.txt").format(sample=sample)),
        targets['assemble'].append(os.path.join(dir_reports, "assembly", "{sample}_assembly_report.txt").format(sample=sample)),
        targets['assemble'].append(os.path.join(dir_reports, "Assembly_stats_all.csv"))


@targetRule
rule all:
    input:
        targets['assemble']