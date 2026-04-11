"""
Mapping rules
For now mapping reads to the individual assemblies
simka rule for picking the 50 assemblies
"""
from glob import glob
sample_names = [s for s in config["sample_names"].keys() if s != "*"]

rule list_simka:
    input:
        r1 = expand(os.path.join(dir_hostcleaned,"{sample}_R1.hostcleaned.fastq.gz"), sample=sample_names),
        r2 = expand(os.path.join(dir_hostcleaned,"{sample}_R2.hostcleaned.fastq.gz"), sample=sample_names)
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
        simka -in {input.simka_list}  -max-reads 0 -abundance-min 2 -max-count 100 -max-merge 16 -max-memory {resources.mem_mb} -nb-cores {threads} -out-tmp {params.dir_out}
        mv {params.dir_out}/mat_abundance_jaccard.csv.gz {output.simka_kmerclust}
        """

rule simka_clusters:
    input:
        r1 = os.path.join(dir_hostcleaned,"{sample}_R1.hostcleaned.fastq.gz"),
        r2 = os.path.join(dir_hostcleaned,"{sample}_R2.hostcleaned.fastq.gz"),
        simka_kmerclust = os.path.join(dir_binning, "mat_abundance_jaccard.csv.gz")
    output:
        simka_clusters_txt = os.path.join(dir_binning, "{sample}_cluster_50", "simka_cluster_50.txt")
    params:
        script = os.path.join(dir_script, "simka_similar_samples_script.sh"), 
        sample = "{sample}",
        simka_clusters = os.path.join(dir_binning, "{sample}_cluster_50")
    conda:
        os.path.join(dir_env, "simka.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        bash {params.script} -s {params.sample} -d {input.simka_kmerclust} -n 49 -o {params.simka_clusters}
        touch {output.simka_clusters_txt}
        """

"""
Backmapping:
Map reads from selected samples to each assembly.
"""
