# magbee

## Snakemake workflow to get MAGs from honey bee metagenomes
Still under development! Stable release will be out as a version \

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
        1. dereplicated_

### Workflow 
Magbee options:

```
magbee --help
Usage: magbee [OPTIONS] COMMAND [ARGS]...

  Assembling pure culture phages from both Illumina and Nanopore sequencing
  technology For more options, run: magbee --help

Options:
  -v, --version  Show the version and exit.
  -h, --help     Show this message and exit.

Commands:
  run          Run magbee
  qc           QC magbee
  assembly     Assembly magbee
  backmapping  Backmapping magbee
  binning      Binning magbee
  config       Copy the system default config file
  citation     Print the citation(s) for this toolz
```

### Workflow example commands 

**QC module**
This command runs the steps for quality control using FastP, and host read removal using Minimap2,  `magbee qc`

    
    Usage: magbee qc [OPTIONS] [SNAKE_ARGS].... 
    QC EXAMPLES 
    magbee qc --input <input directory with reads> --extn fq --pattern_r1 _R1 --pattern_r2 _R2 --host_seq <path to host genomes> --sequencing paired --output <output directory> -k --profile slurm
    

**Assembly module** 
This command runs individual assemblies for each sample using megahit and QUAST for contig reports,  `magbee assembly`

    Usage: magbee assembly [OPTIONS] [SNAKE_ARGS]...

    Assembly magbee

    Assembly EXAMPLES 
    magbee assembly --input <input directory with reads> --extn fq --pattern_r1 <fastq.gz> --pattern_r2 <fastq.gz> --sequencing paired --output <output directory> -k --profile slurm

**Backmapping module**
This command maps the reads to the assembled contigs. There are two modes available here,
- individual, where the reads are mapped to the contigs. If there are more than 100 samples, then simka is run to pick the 50 samples
- concatenate, where the assemblies are pooled together, then all the reads are mapped to the pooled assembly

```
    Usage: magbee backmapping [OPTIONS] [SNAKE_ARGS]...

    Backmapping magbee

    Backmapping EXAMPLES 
    magbee backmapping --input <input directory with reads> --extn fq --pattern_r1 <fastq.gz> --pattern_r2 <fastq.gz> --contigs <input directory with contigs> --sequencing paired --mode concatenate --output <output directory> -k
```     

**Binning module**
This command generates the final set of bins. This is done following the below steps 
- Metabat2 and VAMB binning tools
- dasTools generates a final set of bins between the two tools
- CheckM2 of dastool set for completeness and contamination 

```
    Usage: magbee binning [OPTIONS] [SNAKE_ARGS]...

    Binning magbee

    Backmapping EXAMPLES 
    magbee binning --bam_folder <input directory with bamfiles> --contigs <input directory with contigs> --mode individual --output <output directory> -k
    #output folder here, make sure to define it and add absolute path (dastool will error out if relative path is provided)
```
Note: change --mode individual #for bees, or maybe low richness samples this works best

**Species variation module**
This command does 
- dRep to dereplicate the set
- gtdbtk to generate dRep set


**Species taxa module**
This section annotates the reads to 


**Steps still working on**
The --mode concatenate option not tested from binning module
