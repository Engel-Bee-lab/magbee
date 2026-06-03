"""
Rules to evaluate bins quality running Checkm2 and GTDBtk to assign taxonomy.
"""
from glob import glob

rule bins_checkm2_metabat2:
    input:
        bins_dir = os.path.join(dir_binning, "all_metabat2_bins", "done.txt")
    output:
        checkm2_dir = os.path.join(dir_reports, "checkm2_all", "CheckM2_Metabat2_quality_report.tsv")
    params:
        bins=os.path.join(dir_binning, "all_metabat2_bins"),
        outdir=os.path.join(dir_binning, "checkm2_output_metabat2"),
        database = config["databases"]["checkm2_db"]
    conda:
        os.path.join(dir_env, "checkm2.yaml")
    resources:
        mem_mb =config['resources']['medium']['mem_mb'],
        runtime = config['resources']['medium']['runtime']
    threads: 
        config['resources']['medium']['threads']
    shell:
        """
        mkdir -p {params.outdir}
        export CHECKM2_DB_PATH={params.database}
        if ls {params.bins}/*.fa 1> /dev/null 2>&1; then
            checkm2 predict -i {params.bins} -o {params.outdir} --database_path {params.database}/uniref100.KO.1.dmnd \
                -x fa --force --threads {threads}
            
            if [ -f {params.outdir}/checkm2_assessment.tsv ]; then
                cp {params.outdir}/checkm2_assessment.tsv {output.checkm2_dir}
            else
                echo "CheckM2 failed for metabat2 bins" > {output.checkm2_dir}
            fi
        else
            echo "No bins found, skipping CheckM2 for metabat2 bins" > {output.checkm2_dir}
        fi
        """

rule bins_checkm2_concoct:
    input:
        bins_dir = os.path.join(dir_binning, "all_conoct_bins", "renamed.txt")
    output:
        checkm2_dir = os.path.join(dir_reports, "checkm2_all", "CheckM2_CONCOCT_quality_report.tsv")
    params:
        bins=os.path.join(dir_binning, "all_conoct_bins"),
        outdir=os.path.join(dir_binning, "checkm2_output_concoct"),
        database = config["databases"]["checkm2_db"]
    conda:
        os.path.join(dir_env, "checkm2.yaml")
    resources:
        mem_mb =config['resources']['medium']['mem_mb'],
        runtime = config['resources']['medium']['runtime']
    threads: 
        config['resources']['medium']['threads']
    shell:
        """
        mkdir -p {params.outdir}
        export CHECKM2_DB_PATH={params.database}
        if ls {params.bins}/*.fa 1> /dev/null 2>&1; then
            checkm2 predict -i {params.bins} -o {params.outdir} --database_path {params.database}/uniref100.KO.1.dmnd \
                -x fa --force --threads {threads}
            if [ -f {params.outdir}/checkm2_assessment.tsv ]; then
                cp {params.outdir}/checkm2_assessment.tsv {output.checkm2_dir}
            else
                echo "CheckM2 failed for CONCOCT bins" > {output.checkm2_dir}
            fi
        else
            echo "No bins found, skipping CheckM2 for concoct bins" > {output.checkm2_dir}
        fi
        """

rule bins_checkm2_vamb:
    input:
        bins_dir = os.path.join(dir_binning, "all_vamb_bins", "done.txt")
    output:
        checkm2_dir = os.path.join(dir_reports, "checkm2_all", "CheckM2_VAMB_quality_report.tsv")
    params:
        bins=os.path.join(dir_binning, "all_vamb_bins"),
        outdir=os.path.join(dir_binning, "checkm2_output_vamb"),
        database = config["databases"]["checkm2_db"]
    conda:
        os.path.join(dir_env, "checkm2.yaml")
    resources:
        mem_mb =config['resources']['medium']['mem_mb'],
        runtime = config['resources']['medium']['runtime']
    threads: 
        config['resources']['medium']['threads']
    shell:
        """
        mkdir -p {params.outdir}
        export CHECKM2_DB_PATH={params.database}
        if ls {params.bins}/*.fa 1> /dev/null 2>&1; then
            checkm2 predict -i {params.bins} -o {params.outdir} --database_path {params.database}/uniref100.KO.1.dmnd \
                -x fa --force --threads {threads}
            if [ -f {params.outdir}/checkm2_assessment.tsv ]; then
                cp {params.outdir}/checkm2_assessment.tsv {output.checkm2_dir}
            else
                echo "CheckM2 failed for VAMB bins" > {output.checkm2_dir}
            fi
        else
            echo "No bins found, skipping CheckM2 for vamb bins" > {output.checkm2_dir}
        fi
        """

rule bins_checkm2_semibin:
    input:
        bins_dir = os.path.join(dir_binning, "all_semibin_bins", "done.txt")
    output:
        checkm2_dir = os.path.join(dir_reports, "checkm2_all", "checkM2_SemiBin2_quality_report.tsv")
    params:
        bins=os.path.join(dir_binning, "all_semibin_bins"),
        outdir=os.path.join(dir_binning, "checkm2_output_semibin"),
        database = config["databases"]["checkm2_db"]
    conda:
        os.path.join(dir_env, "checkm2.yaml")
    resources:
        mem_mb =config['resources']['medium']['mem_mb'],
        runtime = config['resources']['medium']['runtime']
    threads: 
        config['resources']['medium']['threads']
    shell:
        """
        mkdir -p {params.outdir}
        export CHECKM2_DB_PATH={params.database}
        if ls {params.bins}/*.fa 1> /dev/null 2>&1; then
            checkm2 predict -i {params.bins} -o {params.outdir} --database_path {params.database}/uniref100.KO.1.dmnd \
                -x fa --force --threads {threads}
            if [ -f {params.outdir}/checkm2_assessment.tsv ]; then
                cp {params.outdir}/checkm2_assessment.tsv {output.checkm2_dir}
            else
                echo "CheckM2 failed for Semibin2 bins" > {output.checkm2_dir}
            fi
        else
            echo "No bins found, skipping CheckM2 for Semibin2 bins" > {output.checkm2_dir}
        fi
        """


