
"""
InStrain rule for Snakemake workflow
"""

#sample_names = [s for s in config["sample_names"].keys() if s != "*"]
sample_names = list(config.get("sample_names", {}).keys())

# Constrain wildcards to actual sample names
wildcard_constraints:
    sample="|".join(sample_names)

rule map_reads2mags:
    input:
        drep_dir = os.path.join(dir_species, "drep_dastool", "done.txt"),
        r1 = os.path.join(dir_hostcleaned,"{sample}_R1.hostcleaned.fastq.gz"),
        r2 = os.path.join(dir_hostcleaned,"{sample}_R2.hostcleaned.fastq.gz")
    params:
        genomes_dir = os.path.join(dir_species, "drep_dastools"),
    output:
        os.path.join(dir_species, "inStrain", "{sample}", "mapped.bam")
    conda:
        os.path.join(dir_env, "minimap2.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    