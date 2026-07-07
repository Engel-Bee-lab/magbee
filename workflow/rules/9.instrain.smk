"""
InStrain workflow for dereplicated MAGs
"""

import glob
import os

sample_names = list(config.get("sample_names", {}).keys())
derep_bins_dir = os.path.join(dir_species, "drep_dastools", "dereplicated_genomes")

def get_derep_bin_files():
    return sorted(glob.glob(os.path.join(derep_bins_dir, "*.fa*")))

def _strip_derep_suffixes(path):
    name = os.path.basename(path)
    for suffix in (".gz", ".fasta", ".fa", ".fna"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name

rule make_mag_rep_database:
    input:
        summary = os.path.join(dir_species, "drep_dastools", "done.txt"),
    output:
        mag_rep_database = os.path.join(dir_species, "inStrain", "prepare_mags", "mag_rep_database.fa"),
    params:
        rep_dir = derep_bins_dir,
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.mag_rep_database})
        cat {params.rep_dir}/*.fa* > {output.mag_rep_database}
        """

rule make_scaffold_to_bin_file:
    input:
        summary = os.path.join(dir_species, "drep_dastools", "done.txt"),
    output:
        scaffold_to_bin_file = os.path.join(dir_species, "inStrain", "prepare_mags", "scaffold_to_bin_file.tsv"),
    params:
        rep_dir = derep_bins_dir,
    run:
        os.makedirs(os.path.dirname(output.scaffold_to_bin_file), exist_ok=True)
        with open(output.scaffold_to_bin_file, "w") as out_f:
            for fasta_path in sorted(glob.glob(os.path.join(params.rep_dir, "*.fa*"))):
                mag_name = _strip_derep_suffixes(fasta_path)
                with open(fasta_path) as fasta_fh:
                    for line in fasta_fh:
                        if line.startswith(">"):
                            scaffold = line[1:].strip().split()[0]
                            out_f.write(f"{scaffold}\t{mag_name}\n")

rule minimp2_magDB_index:
    input:
        mag_rep_database = os.path.join(dir_species, "inStrain", "prepare_mags", "mag_rep_database.fa"),
    output:
        index_done = os.path.join(dir_species, "inStrain", "prepare_mags", "mag_rep_database_minimap2.mmi"),
    conda:
        os.path.join(dir_env, "minimap2.yaml")
    resources:
        mem_mb = config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads:
        config['resources']['smalljob']['threads']
    shell:
        """
        minimap2 -d {output.index_done} {input.mag_rep_database}
        """

rule map_to_rep_MAGs_minimap2:
    input:
        reads1 = lambda wildcards: os.path.join(dir_hostcleaned, f"{wildcards.sample}_R1.hostcleaned.fastq.gz"),
        reads2 = lambda wildcards: os.path.join(dir_hostcleaned, f"{wildcards.sample}_R2.hostcleaned.fastq.gz"),
        mag_rep_database = os.path.join(dir_species, "inStrain", "prepare_mags", "mag_rep_database.fa"),
        index_done = os.path.join(dir_species, "inStrain", "prepare_mags", "mag_rep_database_minimap2.mmi"),
    conda:
        os.path.join(dir_env, "minimap2.yaml")
    output:
        bam = os.path.join(dir_species, "inStrain", "mapping", "{sample}", "{sample}_bowtie.bam"),
        flagstat = os.path.join(dir_species, "inStrain", "mapping", "{sample}", "{sample}_bowtie_flagstat.tsv"),
    params:
        outdir = os.path.join(dir_species, "inStrain", "mapping", "{sample}")
    resources:
        mem_mb = config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads:
        config['resources']['smalljob']['threads']
    shell:
        """
        set -euo pipefail
        if [ -d {params.outdir} ]; then
            rm -rf {params.outdir}
            mkdir -p {params.outdir}
        else
            mkdir -p {params.outdir}
        fi
        minimap2 -ax sr -t {threads} {input.index_done} {input.reads1} {input.reads2} \
            | samtools view -bh - | samtools sort - > {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        """


rule instrain_profile_db_mode:
    input:
        bam = os.path.join(dir_species, "inStrain", "mapping", "{sample}", "{sample}_bowtie.bam"),
        mag_rep_database = os.path.join(dir_species, "inStrain", "prepare_mags", "mag_rep_database.fa"),
        scaffold_to_bin_file = os.path.join(dir_species, "inStrain", "prepare_mags", "scaffold_to_bin_file.tsv"),
    output:
        marker = os.path.join(dir_species, "inStrain", "instrain_profile_db_mode", "{sample}_profile_db_mode.done"),
    params:
        outdir = os.path.join(dir_species, "inStrain", "instrain_profile_db_mode", "{sample}"),
    conda:
        os.path.join(dir_env, "instrain.yaml")
    resources:
        mem_mb = config['resources']['medium']['mem_mb'],
        runtime = config['resources']['medium']['runtime']
    threads:
        config['resources']['medium']['threads']
    shell:
        """
        set -euo pipefail
        module load samtools
        mkdir -p {params.outdir}
        inStrain profile {input.bam} {input.mag_rep_database} -o {params.outdir} \
            -p {threads} --database_mode -s {input.scaffold_to_bin_file}
        touch {output.marker}
        """


if sample_names:

    rule instrain_compare:
        input:
            markers = expand(os.path.join(dir_species, "inStrain", "instrain_profile_db_mode", "{sample}_profile_db_mode.done"), sample=sample_names),
            scaffold_to_bin_file = os.path.join(dir_species, "inStrain", "prepare_mags", "scaffold_to_bin_file.tsv"),
        output:
            compare_marker = touch(os.path.join(dir_species, "inStrain", "instrain_compare", "all_compared.done")),
        params:
            outdir = os.path.join(dir_species, "inStrain", "instrain_compare", "MAGs_rep_db"),
            profiles = [os.path.join(dir_species, "inStrain", "instrain_profile_db_mode", sample) for sample in sample_names],
        conda:
            os.path.join(dir_env, "instrain.yaml")
        resources:
            mem_mb = config['resources']['highmemjob']['mem_mb'],
            runtime = config['resources']['highmemjob']['runtime']
        threads:
            config['resources']['highmemjob']['threads']
        shell:
            """
            set -euo pipefail
            mkdir -p {params.outdir}
            inStrain compare -i {params.profiles} -s {input.scaffold_to_bin_file} \
                -p {threads} -o {params.outdir} --database_mode || touch {output.compare_marker}.singleton
            touch {output.compare_marker}
            """
