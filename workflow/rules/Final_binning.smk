"""
Rule to get the final stats of the binned contigs to reports
"""
rule final_binning_stats:
    input:
        dastoolout=os.path.join(dir_binning, "das_tool", "dastool_DASTool_summary.tsv"),
    output:
        bins=os.path.join(dir_reports, "bins.done"),
    params:
        dasttool=os.path.join(dir_binning, "das_tool"),
        outdir=os.path.join(dir_reports, "dastool_bins")
    localrule: True
    shell:
        """
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
        python {scripts}/filter_bins.py -i {params.checkm2_dir} -o {params.quality_mags} -s {params.outdir} --ext .fa
        
        touch {output.quality}
        """