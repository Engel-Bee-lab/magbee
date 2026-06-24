"""
Rules to evaluate DAS Tool for building the non-redundant set of bins from the output of the binning tools.
"""
from glob import glob

rule contigs2bin_individual:
    input:
        metabat2_bins = os.path.join(dir_binning, "all_metabat2_bins", "done.txt"),
        vamb_bins = os.path.join(dir_binning, "all_vamb_bins", "done.txt")
    params:
        metabat2_bin_folder = os.path.join(dir_binning, "all_metabat2_bins"),
        vamb_bin_folder = os.path.join(dir_binning, "all_vamb_bins"),
        mk = os.path.join(dir_binning, "das_tool", "scaffolds2bin")
    conda:
        os.path.join(dir_env, "dasttool.yaml")
    localrule: True
    output:
        metabat2= os.path.join(dir_binning, "das_tool", "scaffolds2bin", "metabat2_scaffolds2bin.tsv"),
        vamb= os.path.join(dir_binning, "das_tool", "scaffolds2bin", "vamb_scaffolds2bin.tsv")
    shell:
        """
        mkdir -p {params.mk}
        set +e
        Fasta_to_Contig2Bin.sh -i {params.metabat2_bin_folder}/* -e fa > {output.metabat2}

        Fasta_to_Contig2Bin.sh -i {params.vamb_bin_folder}/* -e fasta > {output.vamb}
        """

rule run_DAS_tool_individual:
    input:
        metabat2= os.path.join(dir_binning, "das_tool", "scaffolds2bin", "metabat2_scaffolds2bin.tsv"),
        vamb= os.path.join(dir_binning, "das_tool", "scaffolds2bin", "vamb_scaffolds2bin.tsv"),
        contigs= expand(os.path.join(dir_assembly,"{sample}.megahit.contigs.fa.gz"), sample=sample_names)
    params:
        basename= "dastool",
        sample= "{sample}",
        temp= os.path.join(dir_binning, "das_tool", "temp"),
        bins_dir= os.path.join(dir_binning, "das_tool", "dastool_DASTool_bins")
    conda:
        os.path.join(dir_env, "dasttool.yaml")
    output:
        out=os.path.join(dir_binning, "das_tool", "dastool_DASTool_summary.txt"),
        bins_done=os.path.join(dir_binning, "das_tool", "dastool_DASTool_bins.done")
    threads: 4
    shell:
        """
        zcat {input.contigs} | sed '/^>/s/>/>{params.sample}_/'  > {params.temp}.contigs.fa

        DAS_Tool -i {input.metabat2},{input.vamb} \
            -c {params.temp}.contigs.fa -o {params.basename} --threads {threads} \
            --labels metabat2,vamb  --write_bin_evals --write_bins

        touch {output.bins_done}
        """