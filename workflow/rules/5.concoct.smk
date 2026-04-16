"""
Binning rule using coconct
"""
from glob import glob

rule cut_up_fasta_simka:
    input:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa"),
        bam=os.path.join(dir_binning, "{sample}_cluster_50", "done.txt")
    output:
        bed=os.path.join(dir_binning, "{sample}_concoct", "{sample}.bed"),
        contigs10k=os.path.join(dir_binning, "{sample}_concoct", "{sample}.contigs10k.fa")
    params:
        chunk_size=10000,
        overlap_size=0
    conda:
        os.path.join(dir_env, "concoct.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        cut_up_fasta.py {input.assembly} \
          -c {params.chunk_size} \
          -o {params.overlap_size} \
          --merge_last \
          -b {output.bed} > {output.contigs10K}
        """

rule concoct_table:
    """
    Generates table with per sample coverage depth.
    """
    input:
        bed = os.path.join(dir_binning, "{sample}_concoct", "{sample}.bed"),
        bam=os.path.join(dir_binning, "{sample}_cluster_50", "done.txt")
    output:
        covtable=os.path.join(dir_binning, "{sample}_concoct", "{sample}_covtable.tsv")
    params:
        bam_dir=os.path.join(dir_binning, "{sample}_bam")
    conda:
        os.path.join(dir_env, "concoct.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        concoct_coverage_table.py {input.bed} {params.bam_dir}/*.bam > {output.coverage_table}
        """

rule run_concoct:
    input:
        bed = os.path.join(dir_binning, "{sample}_concoct", "{sample}.bed"),
        contigs10k=os.path.join(dir_binning, "{sample}_concoct", "{sample}.contigs10k.fa"),
        covtable=os.path.join(dir_binning, "{sample}_concoct", "{sample}_covtable.tsv")
    output:
        clusters=os.path.join(dir_binning, "{sample}_concoct", "clustering_gt1000.csv")
    params:
        outdir=os.path.join(dir_binning, "{sample}_concoct"),
        num_clusters=100
    conda:
        os.path.join(dir_env, "concoct.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        concoct --threads {threads} -l 1500 \
            --composition_file {input.contigs10k} --coverage_file {input.covtable} \
            -b {params.outdir} \
            -c {params.num_clusters} \
            -t {threads}
        """

rule concoct_merge:
    input:
        clusters=os.path.join(dir_binning, "{sample}_concoct", "clustering_gt1000.csv")
    output:
        merged=os.path.join(dir_binning, "{sample}_concoct", "clustering_merged.csv")
    conda:
        os.path.join(dir_env, "concoct.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        merge_cutup_clustering.py {input.clusters} > {output.merged}
        """

rule concoct_bins:
    input:
        contigs10k=os.path.join(dir_binning, "{sample}_concoct", "{sample}.contigs10k.fa"),
        merged=os.path.join(dir_binning, "{sample}_concoct", "clustering_merged.csv")
    output:
        bins=os.path.join(dir_binning, "{sample}_concoct", "bins", "done.txt")
    params:
        outdir=os.path.join(dir_binning, "{sample}_concoct", "bins")
    conda:
        os.path.join(dir_env, "concoct.yaml")
    resources:
        mem_mb =config['resources']['long_shortjob']['mem_mb'],
        runtime = config['resources']['long_shortjob']['runtime']
    threads: 
        config['resources']['long_shortjob']['threads']
    shell:
        """
        mkdir -p {params.outdir}
        extract_fasta_bins.py {input.contigs10k} {input.merged} --output_path {params.bins_dir}
        touch {output.bins}
        """

rule rename_concoctBins:
    input:
        bins=os.path.join(dir_binning, "{sample}_concoct", "bins", "done.txt")
    output:
        renamed=os.path.join(dir_binning, "{sample}_concoct", "bins", "renamed.txt")
    params:
        bins_dir=os.path.join(dir_binning, "{sample}_concoct", "bins"),
        sample="{sample}"
    conda:
        os.path.join(dir_env, "concoct.yaml")
    shell:
        """
        for f in {params.bins_dir}/*.fa; do
            base=$(basename $f)
            mv "$f" "{params.bins_dir}/concoct_{params.sample}_$base"
        done
        """