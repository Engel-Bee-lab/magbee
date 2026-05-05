"""
Binning rule using coconct
"""
rule cut_up_fasta_simka:
    input:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa.gz"),
        bam_dir = os.path.join(dir_backmapping, "{sample}_cluster_50", "done.txt")
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
          -b {output.bed} > {output.contigs10k}
        """

rule cut_up_fasta_all2all:
    input:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa.gz"),
        bam_dir = os.path.join(dir_backmapping, "{sample}_bam", "done.txt")
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
          -b {output.bed} > {output.contigs10k}
        """
rule concoct_table_simka:
    """
    Generates table with per sample coverage depth.
    """
    input:
        bed = os.path.join(dir_binning, "{sample}_concoct", "{sample}.bed"),
        bam = os.path.join(dir_backmapping, "{sample}_cluster_50", "done.txt")
    output:
        covtable=os.path.join(dir_binning, "{sample}_concoct", "{sample}_covtable.tsv")
    params:
        bam_dir=os.path.join(dir_binning, "{sample}_cluster_50"),
    conda:
        os.path.join(dir_env, "concoct.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        concoct_coverage_table.py {input.bed} {params.bam_dir}/*.bam > {output.covtable}
        """

rule concoct_table_all2all:
    """
    Generates table with per sample coverage depth.
    """
    input:
        bed = os.path.join(dir_binning, "{sample}_concoct", "{sample}.bed"),
        bam_dir = os.path.join(dir_backmapping, "{sample}_bam", "done.txt")
    output:
        covtable=os.path.join(dir_binning, "{sample}_concoct", "{sample}_covtable.tsv")
    params:
        bam_dir=os.path.join(dir_binning, "{sample}_bam"),
    conda:
        os.path.join(dir_env, "concoct.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        concoct_coverage_table.py {input.bed} {params.bam_dir}/*.bam > {output.covtable}
        """

rule run_concoct:
    input:
        bed = os.path.join(dir_binning, "{sample}_concoct", "{sample}.bed"),
        contigs10k=os.path.join(dir_binning, "{sample}_concoct", "{sample}.contigs10k.fa"),
        covtable=os.path.join(dir_binning, "{sample}_concoct", "{sample}_covtable.tsv")
    output:
        clusters=os.path.join(dir_binning, "{sample}_concoct", "clustering_gt1500.csv")
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
        concoct --threads {threads} -l 1500 --composition_file {input.contigs10k} --coverage_file {input.covtable} \
            -b {params.outdir} -c {params.num_clusters}
        """

rule concoct_merge:
    input:
        clusters=os.path.join(dir_binning, "{sample}_concoct", "clustering_gt1500.csv")
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
        contigs10k=os.path.join(dir_assembly,"{sample}.megahit.contigs.fa.gz"),
        merged=os.path.join(dir_binning, "{sample}_concoct", "clustering_merged.csv")
    output:
        bins=os.path.join(dir_binning, "{sample}_concoct", "bins", "done.txt")
    params:
        outdir=os.path.join(dir_binning, "{sample}_concoct", "bins")
    conda:
        os.path.join(dir_env, "concoct.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        mkdir -p {params.outdir}
        extract_fasta_bins.py {input.contigs10k} {input.merged} --output_path {params.outdir}
        touch {output.bins}
        """

rule rename_concoctBins:
    input:
        bins=expand(os.path.join(dir_binning, "{sample}_concoct", "bins", "done.txt"), sample=sample_names),
    output:
        renamed=os.path.join(dir_binning, "all_conoct_bins", "renamed.txt")
    params:
        outdir=os.path.join(dir_binning, "all_conoct_bins"),
        sample=" ".join(sample_names),
        dir_bins=os.path.join(dir_binning)
    conda:
        os.path.join(dir_env, "concoct.yaml")
    localrule: True
    shell:
        """
        mkdir -p {params.outdir}
        for sample in {params.sample}; do
            src_dir={params.dir_bins}/${{sample}}_concoct/bins
            for f in "$src_dir"/concoct_*.fa; do
                [ -e "$f" ] || continue 
                bn=$(basename "$f")
                cp "$f" "{params.outdir}/${{sample}}_concoct_${{bn}}"
            done
        done
        touch {output.renamed}
        """