"""
Rules to include assembly QUAST output to the reports
"""
from glob import glob

rule assembly_reports:
    input:
        quast_report = os.path.join(dir_assembly, "{sample}_quast_output", "report.txt")
    output:
        report = os.path.join(dir_assembly, "{sample}_assembly_report.txt")
    localrule: True
    shell:
        """
        cp {input.quast_report} {output.report}
        """

rule merge_stats:
    input:
        stats = expand(os.path.join(dir_assembly,"{sample}_assembly_report.txt"), sample=sample_names)
    output:
        merged = os.path.join(dir_reports, "Assembly_stats_all.csv")
    conda:
        os.path.join(dir_env, "scripts.yml")
    params:
        script=os.path.join(dir_script, "merge_quast.py")
    shell:
        """
        python {params.script} -i {input.stats} -o {output.merged}
        """
    
rule one_dir_contigs:
    input:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa")
    output:
        contigs = os.path.join(dir_reports, "assembly", "{sample}.megahit.contigs.fa.gz")
    localrule: True
    params:
        dirs=os.path.join(dir_reports, "assembly")
    shell:
        """
        mkdir -p {params.dirs}
        cp {input.assembly} {params.dirs}/.
        gzip {output.contigs}
        """