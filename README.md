# magbee

## Snakemake workflow to get MAGs from metagenomes
Still under development! Stable release out as a version \
Strain level resolution. 

Documentation: [Read the Docs]() - Working on this 

### Install 

**Source install**
Run the below commands

```
git clone https://github.com/npbhavya/magbee.git
cd magbee
mamba create -y -n magbee python=3.13
conda activate magbee
pip install -e . 
```

conda and pip install will be available when there is stable version release

### Input files
Input files:
Currently supporting only paired-end reads, 
- Input directory with metagenomes
- Reference genome directory
    - include the referene genome for host contamination removal

### Output files
Output files: The user can define the output folder name, or defaults to `magbee.out`; 
- PROCESSING: Folder contains all the intermediate files
- REPORTS: Final results are saved here 
    - QC results: Saved in two parts:
        1. `final_host_qc_summary.txt`: includes a list of metagenomes, total reads, post QC number of reads, number of reads that mapped to referene genome, percent of mapped reads. 
        2. `QC/reads`: This folder contains the reads post QC and host read removal

    - Assembly results: Saved to two parts:
        1. `Assembly_stats_all.csv`: includes a list of the metagenomes and the corresponding quast reports for each metagenome
        2. `assembly`: Asembled contigs.fasta for each metagenome

    - Binning results:

### `magbee run`
The workflow is written to run the whole workflow with `magbee run`
- complete workflow: Runs all the steps
    ```
    magbee run --input metagenome_reads --extn fq.gz --sequencing paired  \
        --pattern_r1 _R1 --pattern_r2 _R2 \
        --host_seq Apis_mellifera_genome_genomic.fna \
        --conda-frontend mamba --output magbee.out \
        --profile slurm
    ```

- host submodule: This part runs fastp and host contamination removal
    ```
    magbee run --input metagenome_reads --extn fq.gz --sequencing paired  \
        --pattern_r1 _R1 --pattern_r2 _R2 \
        --host_seq Apis_mellifera_genome_genomic.fna \
        --conda-frontend mamba --output magbee.out \
        --profile slurm -k host 
    ```

- assembly submodule: QC reads are then used to perform individual assemblies and quast
    ```
    magbee run --input metagenome_reads --extn fq.gz --sequencing paired  \
        --pattern_r1 _R1 --pattern_r2 _R2 \
        --host_seq Apis_mellifera_genome_genomic.fna \
        --conda-frontend mamba --output magbee.out \
        --profile slurm -k assembly 
    ```

- backmapping submodule: This is the part that performs read mapping to the assemblies
There are two strategies here
    1.  