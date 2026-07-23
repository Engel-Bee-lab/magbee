# Example per-sample rule
rule merge_taxa:
    input:
        kraken=os.path.join(dir_taxa, "{sample}", "{sample}.kraken2.report.txt"),
        rat=os.path.join(dir_taxa, "{sample}", "rat.read2classification.abundance.txt")
    output:
        merged=os.path.join(dir_taxa, "{sample}", "{sample}.merged_taxa.tsv")
    conda:
        os.path.join(dir_env, "scripts.yml")
    shell:
        """
        python workflow/scripts/merge_taxa_tables.py \
            --kraken-report {input.kraken} \
            --rat-table {input.rat} \
            --sample {wildcards.sample} \
            --out {output.merged}
        """
        