"""
Rule for running bakta on dereplicated genomes
"""

derep_bins_dir = os.path.join(dir_species, "drep_dastools", "dereplicated_genomes")

rule bakta:
    input:
        summary = os.path.join(dir_species, "drep_dastools", "done.txt"),
    output:
        done = os.path.join(dir_species, "bakta", "done.txt")
    params:
        genome_files = lambda wc: sorted(glob.glob(os.path.join(derep_bins_dir, "*.fa*"))),
        outdir = os.path.join(dir_species, "bakta"),
        db = config["databases"]["bakta_db"]
    conda:
        os.path.join(dir_env, "bakta.yaml")
    resources:
        mem_mb = config['resources']['medium']['mem_mb'],
        runtime = config['resources']['medium']['runtime']
    threads:
        config['resources']['medium']['threads']
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        set -- {params.genome_files}
        if [ "$#" -eq 0 ]; then
            echo "No dereplicated genomes found in {derep_bins_dir}" >&2
            exit 1
        fi

        for genome in "$@"; do
            genome_name=$(basename "$genome")
            genome_name=${{genome_name%.gz}}
            genome_name=${{genome_name%.fa}}
            genome_name=${{genome_name%.fna}}
            genome_name=${{genome_name%.fasta}}

            bakta --db {params.db} --output {params.outdir}/"$genome_name" --threads {threads} \
                "$genome" --prefix "$genome_name" --force
        done

        touch {output.done}
        """
