
"""
InStrain rule for Snakemake workflow
"""

#sample_names = [s for s in config["sample_names"].keys() if s != "*"]
sample_names = list(config.get("sample_names", {}).keys())

# Constrain wildcards to actual sample names
wildcard_constraints:
    sample="|".join(sample_names),
    mag="|".join(mags)

rule map_reads2mags:
    input:
        drep_dir = os.path.join(dir_species, "drep_dastools", "done.txt"),
        genome = lambda wc: next(
            f for f in drep_bin_files
            if extract_mag_name(f) == wc.mag
        ),
        r1 = os.path.join(dir_hostcleaned,"{sample}_R1.hostcleaned.fastq.gz"),
        r2 = os.path.join(dir_hostcleaned,"{sample}_R2.hostcleaned.fastq.gz")
    output:
        os.path.join(dir_species, "inStrain", "{mag}", "{sample}.mapped.bam")
    params:
        os.path.join(dir_species, "inStrain", "{mag}")
    conda:
        os.path.join(dir_env, "minimap2.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        mkdir -p $(dirname {params})
        minimap2 -ax sr -t {threads} {input.genome} {input.r1} {input.r2} | samtools view -bS - > {output}
        """

rule instrain_profile:
    