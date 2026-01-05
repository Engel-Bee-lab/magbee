"""
Rule to get the final stats of the binned contigs to reports
"""
rule final_binning_stats:
    input:
        checkm2_dir = os.path.join(dir_binning, "checkm2_output", "quality_report.tsv"),
        gtdbtk_dir = os.path.join(dir_binning, "gtdbtk_output", "done.txt"),
        drep_dir = os.path.join(dir_binning, "drep_metabat2_bins", "done.txt")
    output:
        final_checkm2_stats = os.path.join(dir_reports, "checkm2_quality_reports.tsv"),
        final_gtdbtk_stats = os.path.join(dir_reports, "gtdbtk_bac120_summary.tsv"),
    params:
        checkm2=os.path.join(dir_binning, "checkm2_output", "quality_report.tsv"),
        gtdbtk=os.path.join(dir_binning, "gtdbtk_output", "gtdbtk.bac120.summary.tsv"),
        drep=os.path.join(dir_binning, "drep_metabat2_bins", "figures"),
        out_figures=os.path.join(dir_reports, "drep_output_figures")
    localrule: True
    shell:
        """
        cp {params.drep} {params.out_figures}
        cp {params.checkm2} {output.final_checkm2_stats}
        cp {params.gtdbtk} {output.final_gtdbtk_stats}
        """