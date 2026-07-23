"""
Taxa rule using Kraken2 to classify host-cleaned paired reads.

This rule consumes the paired host-cleaned FASTQ files validated and
registered by `Reads_taxonomy.smk` via `config["sample_names"]`.
"""

import os

kraken_db = config["databases"].get("kraken_db")
if not kraken_db:
	raise ValueError("Missing config['databases']['kraken_db']")

if not os.path.isdir(kraken_db):
	raise ValueError(f"Kraken2 database not found: {kraken_db}")

rule all:
	input:
		expand(os.path.join(dir_taxa, "{sample}", "kraken2.done"), sample=sample_names)

rule kraken2_reads:
	input:
		r1=lambda wc: config["sample_names"][wc.sample]["r1"],
		r2=lambda wc: config["sample_names"][wc.sample]["r2"],
	output:
		report=os.path.join(dir_taxa, "{sample}", "{sample}.kraken2.report.txt"),
		classification=os.path.join(dir_taxa, "{sample}", "{sample}.kraken2.classifications.tsv"),
		done=os.path.join(dir_taxa, "{sample}", "kraken2.done"),
	params:
		db=kraken_db,
		outdir=lambda wc: os.path.join(dir_taxa, wc.sample),
	threads:
		1
	conda:
		os.path.join(dir_env, "kraken2.yaml")
	shell:
		"""
		set -euo pipefail

		mkdir -p {params.outdir}

		kraken2 \
			--db {params.db} \
			--paired \
			--threads {threads} \
			--report {output.report} \
			--output {output.classification} \
			{input.r1} {input.r2}

		touch {output.done}
		"""