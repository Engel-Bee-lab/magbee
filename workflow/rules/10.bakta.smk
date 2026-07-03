"""
Rule for running bakta on dereplicated genomes
"""

derep_bins_dir = os.path.join(dir_species, "drep_dastools", "dereplicated_genomes")

def get_derep_bin_files():
    return sorted(glob.glob(os.path.join(derep_bins_dir, "*.fa*")))

def _strip_derep_suffixes(filename):
    basename = os.path.basename(filename)
    for suffix in (".gz", ".fasta", ".fa", ".fna"):
        if basename.endswith(suffix):
            basename = basename[: -len(suffix)]
    return basename


def get_derep_bin_names():
    return [_strip_derep_suffixes(path) for path in get_derep_bin_files()]

rule bakta:
    input:
        genome = lambda wc: next(
            path for path in get_derep_bin_files()
            if _strip_derep_suffixes(path) == wc.genome
        ),
    output:
        done = os.path.join(dir_species, "bakta", "{genome}.done")
    params:
        outdir = lambda wc: os.path.join(dir_species, "bakta", wc.genome),
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

        echo "Running bakta on {wildcards.genome}"

        bakta --db {params.db} --output {params.outdir} --prefix {wildcards.genome} --force --threads {threads} \
            {input.genome}

        touch {output.done}
        """
