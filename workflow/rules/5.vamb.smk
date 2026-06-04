"""
Metagenome binning with VAMB
Vamb is a family of metagenomic binners which feeds kmer composition and abundance into a variational autoencoder and clusters the embedding to form bins. 
Its binners perform excellently with multiple samples, and pretty good on single-sample data.
"""
from glob import glob

rule vamb_simka:
    input:
        contigs = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa.gz"),
        bam = os.path.join(dir_backmapping, "{sample}_cluster_50", "done.txt")
    output:
        done = os.path.join(dir_binning, "{sample}_vamb_bins", "vae_clust_done.txt")
    params:
        contigs_rename = os.path.join(dir_binning, "{sample}_vamb_bins", "{sample}_contigs.fa"),
        bin_dir= os.path.join(dir_binning, "{sample}_vamb_bins"),
        bam_dir=os.path.join(dir_backmapping, "{sample}_cluster_50"),
        abundance = os.path.join(dir_binning, "{sample}_vamb_bins", "vae_clusters_unsplit.tsv")
    conda:
        os.path.join(dir_env, "vamb.yaml")
    resources:
        mem_mb =config['resources']['long_shortjob']['mem_mb'],
        runtime = config['resources']['long_shortjob']['runtime']
    threads:
        config['resources']['long_shortjob']['threads']
    shell:
        """
        #using the logic that one assembly was mapped with multiple samples from simka
        if [ -f "{params.abundance}" ]; then
            echo "abundance file already exists, skipping abundance calculation"
            touch {output.done}
        else

            rm -rf {params.bin_dir}
            vamb bin default --outdir {params.bin_dir} --fasta {input.contigs} --bamdir {params.bam_dir} -t {threads}
            touch {output.done}
        fi
        """

rule vamb_all2all:
    input:
        contigs = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa.gz"),
        bam = os.path.join(dir_backmapping, "{sample}_bam", "done.txt")
    output:
        done = os.path.join(dir_binning, "{sample}_vamb_bins", "vae_clust_done.txt")
    params:
        contigs_rename = os.path.join(dir_binning, "{sample}_vamb_bins", "{sample}_contigs.fa"),
        bin_dir= os.path.join(dir_binning, "{sample}_vamb_bins"),
        bam_dir=os.path.join(dir_backmapping, "{sample}_bam"),
        abundance = os.path.join(dir_binning, "{sample}_vamb_bins", "vae_clusters_unsplit.tsv")
    conda:
        os.path.join(dir_env, "vamb.yaml")
    resources:
        mem_mb =config['resources']['long_shortjob']['mem_mb'],
        runtime = config['resources']['long_shortjob']['runtime']
    threads:
        config['resources']['long_shortjob']['threads']
    shell:
        """
        #using the logic that one assembly was mapped with multiple samples from simka
        if [ -f "{params.abundance}" ]; then
            echo "abundance file already exists, skipping abundance calculation"
            touch {output.done}
        else

            rm -rf {params.bin_dir}
            vamb bin default --outdir {params.bin_dir} --fasta {input.contigs} --bamdir {params.bam_dir} -t {threads}
            touch {output.done}
        fi
        """

rule vamb_sep:
    input:
        contigs = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa.gz"),
        bam = os.path.join(dir_binning, "{sample}_vamb_bins", "vae_clust_done.txt")
    output:
        bins = os.path.join(dir_binning, "{sample}_vamb_bins", "done.txt")
    params:
        binsplit=os.path.join(dir_binning, "{sample}_vamb_bins", "vae_clusters_unsplit.tsv"),
        bin_dir= os.path.join(dir_binning, "{sample}_vamb_bins"),
        outdir=os.path.join(dir_binning, "{sample}_vamb_bins", "bins"),
        scripts= os.path.join(dir_script, "vamb_bins_sep.py"),
        min_size=100000
    localrule: True
    shell:
        """
        python {params.scripts} --mapping {params.binsplit} --fasta {input.contigs} \
            --outdir {params.outdir} --min_size {params.min_size}
        rm -rf {params.outdir}/bin_clustername.fasta
        touch {output.bins}
        """

rule vamb_bins:
    input:
        bins = expand(os.path.join(dir_binning, "{sample}_vamb_bins", "done.txt"), sample=sample_names)
    output:
        bins = os.path.join(dir_binning, "all_vamb_bins", "done.txt")
    localrule: True
    params:
        outdir= os.path.join(dir_binning, "all_vamb_bins"),
        sample=" ".join(sample_names),
        dir_binning= os.path.join(dir_binning)
    shell:
        """
        set -euo pipefail

        mkdir -p {params.outdir}

        for sample in {params.sample}; do
            src_dir="{params.dir_binning}/${{sample}}_vamb_bins"

            echo "=== SAMPLE: $sample ==="

            if [ ! -d "$src_dir/bins" ]; then
                echo "Missing bins dir"
                continue
            fi

            n=0

            for f in "$src_dir"/bins/bin_*.fasta; do
                [ -e "$f" ] || continue

                bn=$(basename "$f")
                cp -f "$f" "{params.dir_binning}/all_vamb_bins/${{sample}}_vamb_${{bn}}"

                n=$((n+1))
            done

            echo "Copied $n bins"
        done

        touch {output.bins}
        """

rule bam_dir_make:
    input:
        bam_dirs=os.path.join(dir_backmapping, "{sample}_bam", "done.txt")
    output:
        bam_dir = os.path.join(dir_temp, "all_bam", "{sample}.bam")
    params:
        folder= os.path.join(dir_backmapping, "all_bam"),
        bams=os.path.join(dir_backmapping, "{sample}_bam", "{sample}.bam")
    localrule: True
    shell:
        """
        mkdir -p {params.folder}
        cp {params.bams} {output.bam_dir}
        """

rule vamb_bins_concat:
    input:
        contigs = dir_assembly,
        bam_dir = expand(os.path.join(dir_temp, "all_bam", "{sample}.bam"), sample=sample_names)
    params:
        bin_dir= os.path.join(dir_binning, "vamb_bins_concat"),
        bams=os.path.join(dir_temp, "all_bam"),
        min_size=100000
    output:
        txt = os.path.join(dir_binning, "vamb_bins_concat", "done.txt")
    conda:
        os.path.join(dir_env, "vamb.yaml")
    resources:
        mem_mb = config['resources']['assemblyjob']['mem_mb'],
        runtime = config['resources']['assemblyjob']['runtime']
    threads:
        config['resources']['assemblyjob']['threads']
    shell:
        """
        rm -rf {params.bin_dir}
        vamb bin default --outdir {params.bin_dir} --fasta {input.contigs} --bamdir {params.bams} \
             --minfasta {params.min_size} -o C -m 2000 -t {threads}

        touch {output.txt}
        """