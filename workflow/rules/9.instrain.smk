"""
InStrain workflow for dereplicated MAGs
"""

import glob
import os

def write_wrapped_fasta(source_paths, output_path, line_width=60):
    valid_bases = set("ACGTNacgtn")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w") as out_fh:
        for source_path in source_paths:
            header = None
            sequence_chunks = []

            with open(source_path) as fasta_fh:
                for line_number, raw_line in enumerate(fasta_fh, start=1):
                    line = raw_line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if header is not None:
                            if not sequence_chunks:
                                raise ValueError(
                                    f"{source_path}: header '{header}' has no sequence"
                                )
                            out_fh.write(f">{header}\n")
                            sequence = "".join(sequence_chunks)
                            for start in range(0, len(sequence), line_width):
                                out_fh.write(sequence[start:start + line_width] + "\n")
                        header = line[1:].strip()
                        if not header:
                            raise ValueError(f"{source_path}: empty FASTA header on line {line_number}")
                        sequence_chunks = []
                        continue

                    if header is None:
                        raise ValueError(
                            f"{source_path}: sequence data found before the first header on line {line_number}"
                        )
                    invalid_characters = set(line) - valid_bases
                    if invalid_characters:
                        invalid = "".join(sorted(invalid_characters))
                        raise ValueError(
                            f"{source_path}: invalid character(s) '{invalid}' on line {line_number}"
                        )
                    sequence_chunks.append(line)

            if header is None:
                raise ValueError(f"{source_path}: no FASTA records found")
            if not sequence_chunks:
                raise ValueError(f"{source_path}: header '{header}' has no sequence")
            out_fh.write(f">{header}\n")
            sequence = "".join(sequence_chunks)
            for start in range(0, len(sequence), line_width):
                out_fh.write(sequence[start:start + line_width] + "\n")

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
    run:
        source_paths = sorted(glob.glob(os.path.join(params.rep_dir, "*.fa*")))
        if not source_paths:
            raise ValueError(f"No FASTA files found in {params.rep_dir}")
        write_wrapped_fasta(source_paths, output.mag_rep_database)

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
        #validation_done = os.path.join(dir_species, "inStrain", "prepare_mags", "mag_rep_database.validated"),
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
        bam = os.path.join(dir_species, "inStrain", "mapping", "{sample}", "{sample}_minimap.bam"),
        flagstat = os.path.join(dir_species, "inStrain", "mapping", "{sample}", "{sample}_minimap_flagstat.tsv"),
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
        bam = os.path.join(dir_species, "inStrain", "mapping", "{sample}", "{sample}_minimap.bam"),
        mag_rep_database = os.path.join(dir_species, "inStrain", "prepare_mags", "mag_rep_database.fa"),
        #validation_done = os.path.join(dir_species, "inStrain", "prepare_mags", "mag_rep_database.validated"),
        scaffold_to_bin_file = os.path.join(dir_species, "inStrain", "prepare_mags", "scaffold_to_bin_file.tsv"),
    output:
        marker = os.path.join(dir_species, "inStrain", "instrain_profile_db_mode", "{sample}_profile_db_mode.done"),
    params:
        outdir = os.path.join(dir_species, "inStrain", "instrain_profile_db_mode", "{sample}"),
    conda:
        os.path.join(dir_env, "instrain.yaml")
    resources:
        mem_mb = config['resources']['small_moremem_job']['mem_mb'],
        runtime = config['resources']['small_moremem_job']['runtime']
    threads:
        config['resources']['small_moremem_job']['threads']
    shell:
        """
        set -euo pipefail
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
            mem_mb = config['resources']['small_moremem_job']['mem_mb'],
            runtime = config['resources']['small_moremem_job']['runtime']
        threads:
            config['resources']['small_moremem_job']['threads']
        shell:
            """
            set -euo pipefail
            mkdir -p {params.outdir}
            inStrain compare -i {params.profiles} -s {input.scaffold_to_bin_file} \
                -p {threads} -o {params.outdir} --database_mode || touch {output.compare_marker}.singleton
            touch {output.compare_marker}
            """
