"""
This rule is in progress. 
My idea is to add this in once I have a output I am happy with. 
"""
# Define number of batches (adjust based on your genome count)
NUM_BATCHES = 10 

rule gtdbtk_dastool_individual:
    input:
        bins_dir = os.path.join(dir_binning, "das_tool", "dastool_DASTool_summary.tsv")
    output:
        gtdbtk_dir = os.path.join(dir_binning, "gtdbtk_output_dastool", "done.txt"),
        gtdbtk_summary = os.path.join(dir_binning, "gtdbtk_output_dastool", "classify", "gtdbtk.bac120.summary.tsv"),
        gtdbtk_ar_sumamry = os.path.join(dir_binning, "gtdbtk_output_dastool", "classify", "gtdbtk.ar53.summary.tsv")
    params:
        bins = os.path.join(dir_binning, "das_tool", "dastool_DASTool_bins"),
        outdir = os.path.join(dir_binning, "gtdbtk_dastool_output"),
        database = config["databases"]["gtdbtk_db"]
    conda:
        os.path.join(dir_env, "gtdbtk.yaml")
    resources:
        mem_mb  = config['resources']['assemblyjob']['mem_mb'],
        runtime = config['resources']['assemblyjob']['runtime']
    threads:
        config['resources']['assemblyjob']['threads']
    shell:
        """
        set -euo pipefail

        BATCH_NUM={params.batch_num}
        START=$((BATCH_NUM * 200))
        
        # Create temp batch directory with symlinks
        mkdir -p batch_temp_${{BATCH_NUM}}
        ls {params.bins}/*.fa | tail -n +$((START+1)) | head -n 200 | \
          xargs -I {{}} ln -s $(pwd)/{{}}/.. batch_temp_${{BATCH_NUM}}/
        
        mkdir -p {params.outdir}
        export GTDBTK_DATA_PATH={params.database}
        gtdbtk identify --genome_dir batch_temp_${{BATCH_NUM}} --cpus {threads} --out_dir {params.outdir}/identify -x fa
        gtdbtk align --identify_dir {params.outdir}/identify --out_dir {params.outdir}/align --cpus {threads}
        gtdbtk classify --genome_dir batch_temp_${{BATCH_NUM}} --out_dir {params.outdir}/classify --cpus {threads} -x fa -f --align_dir {params.outdir}/align
        
        # Clean up
        rm -rf batch_temp_${{BATCH_NUM}}
        """

rule combine_gtdbtk_results:
    input:
        bac_summaries = expand(
            os.path.join(dir_binning, "gtdbtk_output_dastool", "batch_{batch_num}", "classify", "gtdbtk.bac120.summary.tsv"),
            batch_num=range(NUM_BATCHES)
        ),
        ar_summaries = expand(
            os.path.join(dir_binning, "gtdbtk_output_dastool", "batch_{batch_num}", "classify", "gtdbtk.ar53.summary.tsv"),
            batch_num=range(NUM_BATCHES)
        )
    output:
        combined_bac = os.path.join(dir_binning, "gtdbtk_output_dastool", "gtdbtk.bac120.summary.tsv"),
        combined_ar = os.path.join(dir_binning, "gtdbtk_output_dastool", "gtdbtk.ar53.summary.tsv"),
        done = os.path.join(dir_binning, "gtdbtk_output_dastool", "done.txt")
    shell:
        """
        # Combine bacterial summaries (skip headers from batches 1+)
        cat {input.bac_summaries[0]} > {output.combined_bac}
        tail -n +2 {input.bac_summaries[1:]} >> {output.combined_bac}
        
        # Combine archaeal summaries
        cat {input.ar_summaries[0]} > {output.combined_ar}
        tail -n +2 {input.ar_summaries[1:]} >> {output.combined_ar}
        
        touch {output.done}
        """

