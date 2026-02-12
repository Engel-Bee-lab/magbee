"""
Mapping rules
For now mapping reads to the individual assemblies
simka rule for picking the 50 assemblies
"""

rule list_simka:
    input:
        r1 = expand(os.path.join(dir_hostcleaned,"{sample}_R1.hostcleaned.fastq.gz"), sample=samples_names),
        r2 = expand(os.path.join(dir_hostcleaned,"{sample}_R2.hostcleaned.fastq.gz"), sample=samples_names)
    output:
        simka_list = os.path.join(dir_binning, "simka_input_list.txt")
    shell:
        """
        # Create the input list for simka
        r1_files=({input.r1})
        r2_files=({input.r2})

        for i in "${!r1_files[@]}"; do
            sample=$(basename "${r1_files[$i]}" | sed 's/_R1.hostcleaned.fastq.gz//')
            echo "${sample}: ${r1_files[$i]} ; ${r2_files[$i]}" >> {output.simka_list}
        done

        """
    
    rule simka_kmerclust:
        input:
            simka_list = os.path.join(dir_binning, "simka_input_list.txt")
        output:
            simka_kmerclust = os.path.join(dir_binning, "mat_abundance_jaccard.csv.gz")
        params:
            dir_out = os.path.join(dir_binning, "simka_temp_output")
        conda:
            os.path.join(dir_env, "simka.yaml")
        resources:
            mem_mb =config['resources']['longjob']['mem_mb'],
            runtime = config['resources']['longjob']['runtime']
        threads: 
            config['resources']['longjob']['threads']
        shell:
            """
            simka -in {input.simka_list}  -max-reads 0 -abundance-min 2 -max-count 100 -max-merge 16 -max-memory {resources.mem_mb} -nb-cores {threads} -out-tmp 
            """