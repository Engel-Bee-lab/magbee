"""
Rules to evaluate DAS Tool for building the non-redundant set of bins from the output of the binning tools.
"""
from glob import glob

rule contigs2bin_individual:
    input:
        metabat2_bins = os.path.join(dir_binning, "all_metabat2_bins", "done.txt"),
        vamb_bins = os.path.join(dir_binning, "all_vamb_bins", "done.txt")
    params:
        metabat2_bin_folder = os.path.join(dir_binning, "all_metabat2_bins"),
        vamb_bin_folder = os.path.join(dir_binning, "all_vamb_bins"),
        mk = os.path.join(dir_binning, "das_tool", "scaffolds2bin")
    conda:
        os.path.join(dir_env, "dasttool.yaml")
    localrule: True
    output:
        metabat2= os.path.join(dir_binning, "das_tool", "scaffolds2bin", "metabat2_scaffolds2bin.tsv"),
        vamb= os.path.join(dir_binning, "das_tool", "scaffolds2bin", "vamb_scaffolds2bin.tsv")
    shell:
        """
        mkdir -p {params.mk}
        set +e
        Fasta_to_Contig2Bin.sh -i {params.metabat2_bin_folder}/* -e fa > {output.metabat2}

        Fasta_to_Contig2Bin.sh -i {params.vamb_bin_folder}/* -e fasta > {output.vamb}
        
        # Post-process to add sample names from bin filenames to first column
        python3 << 'EOF'
import os
import re

def process_scaffolds2bin(tsv_file):
    # Process scaffolds2bin file to keep sample name prefix in bin names
    if not os.path.exists(tsv_file):
        return
    
    with open(tsv_file, 'r') as f:
        lines = f.readlines()
    
    processed_lines = []
    for line in lines:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            contig = parts[0]
            bin_name = parts[1]
            
            # Extract sample name from bin name
            # Assuming bin names are like: sample_bin.1, sample_bin.2, etc.
            match = re.match(r'([^_]+)_(.+)', bin_name)
            if match:
                sample_name = match.group(1)
                bin_id = match.group(2)
                new_bin_name = f"{sample_name}_{bin_id}"
            else:
                new_bin_name = bin_name
            
            processed_lines.append(f"{contig}\t{new_bin_name}\n")
        else:
            processed_lines.append(line)
    
    with open(tsv_file, 'w') as f:
        f.writelines(processed_lines)

process_scaffolds2bin("{output.metabat2}")
process_scaffolds2bin("{output.vamb}")
EOF
        """

rule run_DAS_tool_individual:
    input:
        metabat2= os.path.join(dir_binning, "das_tool", "scaffolds2bin", "metabat2_scaffolds2bin.tsv"),
        vamb= os.path.join(dir_binning, "das_tool", "scaffolds2bin", "vamb_scaffolds2bin.tsv"),
    params:
        basename= "dastool",
        temp_contigs = os.path.join(dir_binning, "das_tool", "temp", "combined.contigs.fa"),
        outdir = os.path.join(dir_binning, "das_tool"),
        bins_dir= os.path.join(dir_binning, "das_tool", "dastool_DASTool_bins")
    conda:
        os.path.join(dir_env, "dasttool.yaml")
    output:
        out=os.path.join(dir_binning, "das_tool", "dastool_DASTool_summary.txt"),
        bins_done=os.path.join(dir_binning, "das_tool", "dastool_DASTool_bins.done")
    threads: 4
    shell:
        """
        mkdir -p {params.outdir}/temp
        
        # Combine contigs from all samples, prefixing each with its sample name
        for contig_file in {input.contigs}; do
            sample_name=$(basename "$contig_file" .megahit.contigs.fa.gz)
            zcat "$contig_file" | sed "/^>/s/>/>${{sample_name}}_/" >> {params.temp_contigs}
        done

        cd {params.outdir}
        DAS_Tool -i {input.metabat2},{input.vamb} \
            -c {params.temp_contigs} -o {params.basename} --threads {threads} \
            --labels metabat2,vamb  --write_bin_evals --write_bins

        touch {output.bins_done}
        """