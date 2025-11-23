"""
Binning rule using metabat2
"""
from glob import glob

rule metabat2_binning_individual_sample:
    input:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa"),
        bam = os.path.join(dir_binning, "{sample}_bam", "done.txt"),
    output:
        bins_dir = os.path.join(dir_binning, "{sample}_metabat2_bins", "done.txt")
    params:
        outdir=os.path.join(dir_binning, "{sample}_metabat2_bins"),
        bam_dir=os.path.join(dir_binning, "{sample}_bam")
    conda:
        os.path.join(dir_env, "metabat2.yaml")
    resources:
        mem =config['resources']['bigjob']['mem'],
        time = config['resources']['bigjob']['time']
    threads: 
        config['resources']['bigjob']['cpu']
    shell:
        """
        mkdir -p {params.outdir}
        jgi_summarize_bam_contig_depths --outputDepth {params.outdir}/depth.txt {params.bam_dir}/*.bam
        metabat2 -i {input.assembly} -a {params.outdir}/depth.txt -o {params.outdir}/bin -t {threads}
        touch {output.bins_dir}
        """

rule drep_metabat2_bins:
    input:
        bins_done = expand(
            os.path.join(dir_binning, "{sample}_metabat2_bins", "done.txt"),
            sample=sample_names
        )
    output:
        drep_dir = os.path.join(dir_binning, "drep_metabat2_bins", "done.txt")
    params:
        outdir      = os.path.join(dir_binning, "drep_metabat2_bins"),
        genomes_dir = os.path.join(dir_binning, "drep_metabat2_bins", "genomes"),
        bins_dirs   = expand(
            os.path.join(dir_binning, "{sample}_metabat2_bins"),
            sample=sample_names
        )
    conda:
        os.path.join(dir_env, "drep.yaml")
    resources:
        mem  = config['resources']['bigjob']['mem'],
        time = config['resources']['bigjob']['time']
    threads:
        config['resources']['bigjob']['cpu']
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.outdir}
        mkdir -p {params.genomes_dir}

        # copy bins with unique names: <sample>_<bin>.fa
        for d in {params.bins_dirs}; do
            sample=$(basename "$d" | sed 's/_metabat2_bins//')
            for f in "$d"/bin*.fa; do
                [ -e "$f" ] || continue
                bn=$(basename "$f")
                newname="${{sample}}_${{bn}}"
                # copy instead of symlink to avoid broken link issues
                cp -f "$f" "{params.genomes_dir}/$newname"
            done
        done

        # sanity check: make sure we actually have genomes
        if ! ls {params.genomes_dir}/*.fa >/dev/null 2>&1; then
            echo "No genome FASTA files found in {params.genomes_dir}" >&2
            exit 1
        fi

        # run dRep on the copied genomes
        dRep dereplicate {params.outdir} \
            -g {params.genomes_dir}/*.fa \
            -p {threads}

        touch {output.drep_dir}
        """
