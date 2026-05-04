mport yaml
import os
import glob

"""Parse config"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")

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