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
        executable("/bin/bash")
        r"""
        set -euo pipefail

        mkdir -p magebee-binning/PROCESSING/5_binning/all_vamb_bins

        for sample in {params.sample}; do
            src_dir="magebee-binning/PROCESSING/5_binning/${{sample}}_vamb_bins"

            echo "=== SAMPLE: $sample ==="

            if [ ! -d "$src_dir/bins" ]; then
                echo "Missing bins dir"
                continue
            fi

            n=0

            for f in "$src_dir"/bins/bin_*.fasta; do
                [ -e "$f" ] || continue

                bn=$(basename "$f")
                cp -f "$f" "magebee-binning/PROCESSING/5_binning/all_vamb_bins/${sample}_vamb_${bn}"

                n=$((n+1))
            done

            echo "Copied $n bins"
        done

        touch magebee-binning/PROCESSING/5_binning/all_vamb_bins/done.txt
        """