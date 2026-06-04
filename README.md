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
**run option**
This command runs all the steps in the whole workflow with `magbee run`.

The options for this command:

    
    magbee run --help
    Usage: magbee run [OPTIONS] [SNAKE_ARGS]...

    Run magbee

    Options:
    --sequencing [paired|longread]  sequencing method  [default: paired]
    --input PATH                    Directory of reads  [default:
                                    testReads/paired]
    --extn PATH                     Reads extension; fastq, fq, fastq.gz
                                    [default: fastq]
    --host_seq PATH                 Path to host genome index for host read
                                    removal
    --pattern_r1 TEXT               Pattern to identify R1 reads (for paired-end
                                    data)  [default: _R1]
    --pattern_r2 TEXT               Pattern to identify R2 reads (for paired-end
                                    data)  [default: _R2]
    --contigs PATH                  Assembled contigs
    --mode [individual|concatenate]
                                    Backmapping mode: "individual" or
                                    "concatenate"  [default: concatenate]
    --bam_folder PATH               Directory of bam files for binning
    --output PATH                   Output directory  [default: output]
    --configfile TEXT               Custom config file [default: config.yaml]
    --threads INTEGER               Number of threads to use  [default: 1]
    --profile TEXT                  Snakemake profile
    --db_dir PATH                   Custom database directory
    --temp-dir TEXT                 Temp directory
    --snake-default TEXT            Customise Snakemake runtime args  [default:
                                    --rerun-incomplete, --printshellcmds,
                                    --nolock, --show-failed-logs]
    --use-conda BOOLEAN             Use conda for Snakemake rules  [default:
                                    True]
    --conda-frontend TEXT           Use mamba for Snakemake rules  [default:
                                    mamba]
    --conda-prefix PATH             Custom conda env directory  [default:
                                    /scratch/bnalaga1/magbee/workflow/conda]
    -h, --help                      Show this message and exit.

    RUN EXAMPLES 
    magbee run --input <input directory with reads> --extn fq --pattern_r1 <fastq.gz> --pattern_r2 <fastq.gz> --host_seq <path to host genomes> --sequencing paired --output <output directory> -k
    

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
This command runs Metabat2 and VAMB binning tools.


    Usage: magbee binning [OPTIONS] [SNAKE_ARGS]...

    Binning magbee

    Backmapping EXAMPLES 
    magbee binning --bam_folder <input directory with bamfiles> --contigs <input directory with contigs> --mode concatenate --output <output directory> -k


#### Which mode works?
