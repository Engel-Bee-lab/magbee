"""
Rule to get the final stats of the binned contigs to reports
"""
rule final_binning_stats:
    input:
        dastoolout=os.path.join(dir_binning, "das_tool", "dastool_DASTool_summary.tsv"),
        gtdbtk_summary = os.path.join(dir_binning, "gtdbtk_output_dastool", "classify", "gtdbtk.bac120.summary.tsv"),
        gtdbtk_ar_sumamry = os.path.join(dir_binning, "gtdbtk_output_dastool", "classify", "gtdbtk.ar53.summary.tsv")
    output:
        final_gtdbtk_bac = os.path.join(dir_reports, "gtdbtk_dastools_bac120_summary.tsv"),
        final_gtdbtk_arc = os.path.join(dir_reports, "gtdbtk_dastools_ar53_summary.tsv"),
        bins=os.path.join(dir_reports, "bins.done"),
    params:
        gtdbtk=os.path.join(dir_binning, "gtdbtk_output", "classify"),
        dasttool=os.path.join(dir_binning, "das_tool"),
        outdir=os.path.join(dir_reports, "dastool_bins")
    localrule: True
    shell:
        """
        cp {input.gtdbtk_summary} {output.final_gtdbtk_bac}
        cp {input.gtdbtk_ar_summary} {output.final_gtdbtk_arc}
        mv {params.dasttool}/dastool_DASTool_bins {params.outdir}
        touch {output.bins}
        """
    
rule quality_mags:
    input:
        dasttool=os.path.join(dir_reports, "bins.done"),
        checkm2_dir = os.path.join(dir_reports, "CheckM2_DASTool_quality_report.tsv")
    output:
        quality=os.path.join(dir_reports, "quality_mags.done"),
    params:
        outdir=os.path.join(dir_reports, "dastool_bins"),
        quality_mags=os.path.join(dir_reports, "quality_mags"),
        scripts= os.path.join(dir_script, "filter_bins.py"),
    localrule: True
    shell:
        """
        python {scripts}/filter_bins.py -i {params.checkm2_dir} -o {params.quality_mags} \
            --hq_comp 95 --hq_cont 5 \
            --mq_comp_min 50 --mq_comp_max 90 --mq_cont 10
        
        touch {output.quality}
        """