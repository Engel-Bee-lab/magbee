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
        # Generate TSV and add sample prefix to contig names
        Fasta_to_Contig2Bin.sh -i {params.metabat2_bin_folder} -e fa > {output.metabat2}

        Fasta_to_Contig2Bin.sh -i {params.vamb_bin_folder} -e fasta  > {output.vamb}
        """

rule run_DAS_tool_individual:
    input:
        metabat2= os.path.join(dir_binning, "das_tool", "scaffolds2bin", "metabat2_scaffolds2bin.tsv"),
        vamb= os.path.join(dir_binning, "das_tool", "scaffolds2bin", "vamb_scaffolds2bin.tsv"),
        contigs= expand(os.path.join(dir_assembly,"{sample}.annotated.megahit.contigs.fa.gz"), sample=sample_names)
    params:
        basename= "dastool",
        temp_contigs = os.path.join(dir_binning, "das_tool", "temp", "combined.contigs.fa"),
        outdir = os.path.join(dir_binning, "das_tool"),
        bins_dir= os.path.join(dir_binning, "das_tool", "dastool_DASTool_bins")
    conda:
        os.path.join(dir_env, "dasttool.yaml")
    output:
        out=os.path.join(dir_binning, "das_tool", "dastool_DASTool_summary.tsv"),
    threads: 4
    shell:
        """
        mkdir -p {params.outdir}/temp
        
        zcat {input.contigs} >> {params.temp_contigs}

        cd {params.outdir}
        DAS_Tool -i {input.metabat2},{input.vamb} \
            -c {params.temp_contigs} -o {params.basename} --threads {threads} \
            --labels metabat2,vamb  --write_bin_evals --write_bins
        """