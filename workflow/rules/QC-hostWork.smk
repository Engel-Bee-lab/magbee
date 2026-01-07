"""
Extra steps with the host reads 
    - map to the the reference mitogenome
"""
from glob import glob

rule host_mito_mapping:
    input:
        mr1 = os.path.join(dir_hostcleaned,"{sample}_R1.mapped.fastq.gz"),
        mr2 = os.path.join(dir_hostcleaned,"{sample}_R2.mapped.fastq.gz"),
        host= config['extra_db']['mitogenome']
    output:
        mr1_mt = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mt_R1.mapped.fastq.gz"),
        mr2_mt = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mt_R2.mapped.fastq.gz"),
    params:
        mapped_bam=os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mapped.bam"),
        stats=os.path.join(dir_hostcleaned, "mitogenome", "{sample}_bamstats.txt")
    conda:
        os.path.join(dir_env, "minimap2.yaml")
    resources:
        mem =config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    threads: 
        config['resources']['smalljob']['cpu']
    shell:
        """
        set -euo pipefail
        if [ -f {output.mr1_mt} ] && [ -f {output.mr2_mt} ]; then
            echo "Host mitogenome mapping already done. Skipping..."
            exit 0
        else
            minimap2 -ax sr -t {threads} {input.host} {input.mr1} {input.mr2} \
                | samtools view -b -F 4 -@ {threads} -o {params.mapped_bam} -
            samtools flagstat {params.mapped_bam} > {params.stats}

            samtools fastq -@ {threads} -0 /dev/null -s /dev/null -n \
                -1 >(gzip -c  > {output.mr1_mt}) \
                -2 >(gzip -c  > {output.mr2_mt}) \
                {params.mapped_bam}
            

            touch {output.mr1_mt}
            touch {output.mr2_mt}
        
        fi
        """

rule host_mito_snps:
    input:
        bam = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mapped.bam"),
        host= config['extra_db']['mitogenome']
    output:
        vcf = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps.vcf.gz"),
        filterred_vcf = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps.filtered.vcf.gz")
    params:
        sort_bam=os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mapped.sorted.bam"),
        prefix = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome"),
        stats=os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps.stats.txt"),
        depth=os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps_depth.txt")
    conda:
        os.path.join(dir_env, "bcftools.yaml")
    resources:
        mem =config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    threads: 
        config['resources']['smalljob']['cpu']
    shell:
        """
        set -euo pipefail
        if [ -f {output.vcf} ]; then
            echo "Output vcf file found, so this run looks liks its run. Skipping..."
            exit 0
        else

            #prep the bam
            samtools sort -o {params.sort_bam} {input.bam}
            samtools index {params.sort_bam}

            #call variants
            bcftools mpileup -f {input.host} -Q 20 -q 20 {params.sort_bam} | \
                bcftools call --ploidy 1 -mv -Ov -o {output.vcf}
            bcftools index {output.vcf}

            #filter the SNPs, QUAL>30 means 0.1% error rate, DP>10 means at least 10 reads support
            bcftools filter -i 'QUAL>30 && DP>10' {output.vcf} -Oz -o {output.filterred_vcf}
            bcftools index {output.filterred_vcf}

            #getting the stats
            bcftools stats {output.filterred_vcf} > {params.stats}
            samtools depth {params.sort_bam} > {params.depth}
        fi
        """

rule normalise_vcfs:
    input:
        filterred_vcf = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps.filtered.vcf.gz")
    output:
        norm_vcf = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps.filtered.norm.vcf.gz")
    params:
        host= config['extra_db']['mitogenome'],
        temp=os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps.filtered.temp.vcf.gz"), 
    conda:
        os.path.join(dir_env, "bcftools.yaml")
    shell:
        """
        set -euo pipefail
        if [ -f {output.norm_vcf} ]; then
            echo "Normalized VCF already exists. Skipping..."
            exit 0
        else
            bcftools norm -f {params.host} -m -any {input.filterred_vcf} -Oz -o {params.temp} 
            bcftools index {params.temp}

            #then keep only SNPs
            bcftools view -v snps {params.temp} -Oz -o {output.norm_vcf}
            bcftools index {output.norm_vcf}
        fi
        """

rule merge_vcf:
    input:
        expand(os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps.filtered.norm.vcf.gz"), sample=sample_names)
    output:
        merged_vcf = os.path.join(dir_hostcleaned, "mitogenome", "merged_mitogenome_snps.filtered.norm.vcf.gz")
    params:
        temp=os.path.join(dir_hostcleaned, "mitogenome", "merged_mitogenome_snps.temp.vcf.gz"),
        folder=os.path.join(dir_hostcleaned, "mitogenome")
    conda:
        os.path.join(dir_env, "bcftools.yaml")
    shell:
        """
        set -euo pipefail
        if [ -f {output.merged_vcf} ]; then
            echo "Merged VCF already exists. Skipping..."
            exit 0
        else
            bcftools merge -Oz -o {params.temp} {params.folder}/*_mitogenome_snps.filtered.norm.vcf.gz

            mv {params.temp} {output.merged_vcf}
            bcftools index {output.merged_vcf}
        fi
        """

rule snp_alignment:
    input:
        vcf = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps.filtered.norm.vcf.gz"),
        merged_vcf = os.path.join(dir_hostcleaned, "mitogenome", "merged_mitogenome_snps.filtered.norm.vcf.gz"),
        host= config['extra_db']['mitogenome']
    output:
        consensus_fasta = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_consensus.fasta")
    params:
        sample = "{sample}",
    conda:
        os.path.join(dir_env, "bcftools.yaml")
    params:
        sample="{sample}"
    shell:
        """
        set -euo pipefail
        if [ -f {output.consensus_fasta} ]; then
            echo "SNP alignment fasta already exists. Skipping..."
            exit 0
        else
            bcftools consensus -s {wildcards.sample} -f {input.host} {input.merged_vcf} > {output.consensus_fasta}
        fi
        """

rule build_alignment_fasta:
    input:
        fasta=expand(os.path.join(dir_hostcleaned, "mitogenome", "{sample}_consensus.fasta"), sample=sample_names)
    output:
        final_fasta = os.path.join(dir_hostcleaned, "mitogenome", "final_mitogenome_alignment.fasta")
    params:
        folder=os.path.join(dir_hostcleaned, "mitogenome")
    conda:
        os.path.join(dir_env, "bcftools.yaml")
    shell:
        """
        set -euo pipefail
        if [ -f {output.final_fasta} ]; then
            echo "Final alignment fasta already exists. Skipping..."
            exit 0
        else
            snp-sites -o {output.final_fasta} {params.folder}/*_consensus.fasta
        fi
        """

rule phylo_tree:
    input:
        final_fasta = os.path.join(dir_hostcleaned, "mitogenome", "final_mitogenome_alignment.fasta")
    output:
        tree = os.path.join(dir_reports, "mitogenome", "mitogenome_phylo_tree.nwk")
    conda:
        os.path.join(dir_env, "bcftools.yaml")
    shell:
        """
        set -euo pipefail
        if [ -f {output.tree} ]; then
            echo "Phylogenetic tree already exists. Skipping..."
            exit 0
        else
            iqtree -s {input.final_fasta} -m GTR+G -bb 3000 -nt AUTO
        fi
        """