"""
Rule to get the final backmapping outputs 
"""

rule report_all_to_all:
    input:
        samples= os.path.join(dir_backmapping, "{sample}_bam", "done.txt")
    output:
        os.path.join(dir_reports, "backmapping_report_all_to_all.txt")
    params:
        folder=dir_backmapping
    localrule:True
    shell:
        """
        echo "All-to-All mapping strategy" >> {output}
        echo "Mapped bam files saved to {params.folder}" >> {output}
        echo "These files take a lot of space, so they are not saved in the final output folder." >> {output}
        echo "Informing where these files can be found for the next module, if you are running this workflow modularly." >> {output}
        """

rule report_simka_strategy:
    input:
        done = os.path.join(dir_backmapping, "{sample}_cluster_50", "done.txt")
    output:
        os.path.join(dir_reports, "backmapping_report_simka.txt")
    params:
        folder=dir_backmapping
    localrule:True
    shell:
        """
        echo "Simka mapping strategy" >> {output}
        echo "Mapped bam files saved to {params.folder}" >> {output}
        echo "These files take a lot of space, so they are not saved in the final output folder." >> {output}
        echo "Informing where these files can be found for the next module, if you are running this workflow modularly." >> {output}
        """