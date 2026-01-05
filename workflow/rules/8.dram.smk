rule dram:
    input:
        bins_done = os.path.join(dir_binning, "all_metabat2_bins", "done.txt")
    output:
        dram_out = directory("annotation/dram_annotation")
    params:
        dram_db = config["dram_db"]
    threads: 8
    shell:
        """
        dram.py annotate -i {input.dram_fa} -g {input.dram_gtf} \
            -o {output.dram_out} --threads {threads} --db_dir {params.dram_db}
        """