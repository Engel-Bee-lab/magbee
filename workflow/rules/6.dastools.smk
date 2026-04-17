"""
Rules to evaluate DAS Tool for building the non-redundant set of bins from the output of the binning tools.
"""

rule run_DAS_tool:
    input:
        
    output:
        
    threads: 4
    shell:
        """
        DAS_Tool -i {input.binning_output} -c {input.contigs} -o {output.das_tool_output} --threads {threads} {params.das_tool_params}
        """