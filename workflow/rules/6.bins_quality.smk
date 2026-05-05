"""
Rules to evaluate bins quality running Checkm2 and GTDBtk to assign taxonomy.
"""
from glob import glob

rule bins_checkm2_metabat2:
    input:
        bins_dir = os.path.join(dir_binning, "all_metabat2_bins", "done.txt")
    output:
        checkm2_dir = os.path.join(dir_binning, "checkm2", "checkm2_output_metabat2", "quality_report.tsv")
    params:
        bins=os.path.join(dir_binning, "all_metabat2_bins"),
        outdir=os.path.join(dir_binning, "checkm2_output_metabat2"),
        database = config["databases"]["checkm2_db"]
    conda:
        os.path.join(dir_env, "checkm2.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        mkdir -p {params.outdir}
        export CHECKM2_DB_PATH={params.database}
        if ls {params.bins}/*.fa 1> /dev/null 2>&1; then
            checkm2 predict -i {params.bins} -o {params.outdir} --database_path {params.database}/uniref100.KO.1.dmnd \
                -x fa --force --threads {threads}
            cp {params.outdir}/checkm2_assessment.tsv {output.checkm2_dir}
        else
            echo "No bins found, skipping CheckM2 for metabat2 bins" > {output.checkm2_dir}
        fi
        """

rule bins_checkm2_concoct:
    input:
        bins_dir = os.path.join(dir_binning, "all_conoct_bins", "renamed.txt")
    output:
        checkm2_dir = os.path.join(dir_binning, "checkm2", "checkm2_output_concoct", "quality_report.tsv")
    params:
        bins=os.path.join(dir_binning, "all_conoct_bins"),
        outdir=os.path.join(dir_binning, "checkm2_output_concoct"),
        database = config["databases"]["checkm2_db"]
    conda:
        os.path.join(dir_env, "checkm2.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        mkdir -p {params.outdir}
        export CHECKM2_DB_PATH={params.database}
        if ls {params.bins}/*.fa 1> /dev/null 2>&1; then
            checkm2 predict -i {params.bins} -o {params.outdir} --database_path {params.database}/uniref100.KO.1.dmnd \
                -x fa --force --threads {threads}
            cp {params.outdir}/checkm2_assessment.tsv {output.checkm2_dir}
        else
            echo "No bins found, skipping CheckM2 for concoct bins" > {output.checkm2_dir}
        fi
        """

rule bins_checkm2_vamb:
    input:
        bins_dir = os.path.join(dir_binning, "all_vamb_bins", "done.txt")
    output:
        checkm2_dir = os.path.join(dir_binning, "checkm2", "checkm2_output_vamb", "quality_report.tsv")
    params:
        bins=os.path.join(dir_binning, "all_vamb_bins"),
        outdir=os.path.join(dir_binning, "checkm2_output_vamb"),
        database = config["databases"]["checkm2_db"]
    conda:
        os.path.join(dir_env, "checkm2.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        mkdir -p {params.outdir}
        export CHECKM2_DB_PATH={params.database}
        if ls {params.bins}/*.fa 1> /dev/null 2>&1; then
            checkm2 predict -i {params.bins} -o {params.outdir} --database_path {params.database}/uniref100.KO.1.dmnd \
                -x fa --force --threads {threads}
            cp {params.outdir}/checkm2_assessment.tsv {output.checkm2_dir}
        else
            echo "No bins found, skipping CheckM2 for vamb bins" > {output.checkm2_dir}
        fi
        """

rule gtdbtk_bins:
    input:
        bins_done = os.path.join(dir_binning, "all_metabat2_bins", "done.txt"),
    output:
        gtdbtk_dir = os.path.join(dir_binning, "gtdbtk_output", "done.txt")
    params:
        bins = os.path.join(dir_binning, "all_metabat2_bins"),
        outdir = os.path.join(dir_binning, "gtdbtk_output"),
        database = config["databases"]["gtdbtk_db"]
    conda:
        os.path.join(dir_env, "gtdbtk.yaml")
    resources:
        mem_mb  = config['resources']['bigjob']['mem_mb'],
        runtime = config['resources']['bigjob']['runtime']
    threads:
        config['resources']['bigjob']['threads']
    shell:
        """
        set -euo pipefail

        mkdir -p {params.outdir}
        export GTDBTK_DATA_PATH={params.database}
        gtdbtk identify --genome_dir {params.bins} --cpus {threads} --out_dir {params.outdir} -x fa
        gtdbtk align --identify_dir {params.outdir} --out_dir {params.outdir} --cpus {threads} 
        gtdbtk classify --genome_dir {params.bins} --out_dir {params.outdir} --cpus {threads} -x fa -f --align_dir {params.outdir}
        touch {output.gtdbtk_dir}
        """



