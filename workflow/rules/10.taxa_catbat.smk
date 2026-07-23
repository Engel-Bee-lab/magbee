"""
CATBATRAT workflow
The idea here is to get the the drep bin info to annotate the reads using RAT workflow
But also maybe run kraken2 to get a better idea of the taxonomy of the reads
"""

import glob
import os

rat_db_root = config['databases'].get('rat_db_dir', None)
if not rat_db_root:
    raise ValueError("Missing rat_db_dir in config. Provide the CAT_pack database root containing db/ and tax/.")

rat_db_dir = os.path.join(rat_db_root, "db")
rat_tax_dir = os.path.join(rat_db_root, "tax")

rule rat:
    input:
        r1 = lambda wc: config["sample_names"][wc.sample]["r1"],
        r2 = lambda wc: config["sample_names"][wc.sample]["r2"],
        contigs = config["args"]["contigs"],
        bins = config["args"]["bins"],
    output:
        done = os.path.join(dir_taxa, "{sample}", "done.txt"),
    params:
        out_prefix = lambda wc: os.path.join(dir_taxa, wc.sample, "rat"),
        db_dir = rat_db_dir,
        tax_dir = rat_tax_dir,
    threads:
        1
    conda:
        os.path.join(dir_env, "catbatrat.yaml")
    shell:
        """
        set -euo pipefail

        mkdir -p {params.out_prefix}

        CAT_pack reads --mode mcr \
            -b {input.bins} \
            -c {input.contigs} \
            -1 {input.r1} \
            -2 {input.r2} \
            -d {params.db_dir} \
            -t {params.tax_dir} \
            -o {params.out_prefix}

        touch {output.done}
        """
