#Snakemake rule to amke sure the input files are correct and all the params are provided to start the workflow

#To start with the input files have to be in a folders
#can be compressed or decompressed
#For now there can be no sub folders, the files have to be saved in the same folder

import os 
import glob
import yaml
import re
from metasnek import fastq_finder

"""
CONFIG FILE
"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")

"""
DIRECTORIES
"""
dir = {}
#declaring output file
try:
    if config['args']['output'] is None:
        dir_out = os.path.join('output')
    else:
	    dir_out = config['args']['output']
except KeyError:
    dir_out = os.path.join('output')

# temp dir
if config['args']['temp_dir'] is None:
    dir_temp = os.path.join(dir_out, "temp")
else:
    dir_temp = config['args']['temp_dir']


#declaring some the base directories
dir_env = os.path.join(workflow.basedir,"envs")
dir_script = os.path.join(workflow.basedir,"scripts")

#making directories for each step
#Saving most of the files to PROCESSING, sine they are intermediate files
dir_fastp = os.path.join(dir_out, 'PROCESSING' ,'1_fastp')
dir_hostcleaned = os.path.join(dir_out, 'PROCESSING' ,'2_host_cleaned')
dir_reports = os.path.join(dir_out, 'REPORTS')
#dir_assembly = os.path.join(dir_out, 'PROCESSING' ,'3_coassembly')
dir_assembly = os.path.join(dir_out, 'PROCESSING' ,'3_individual_assembly')
dir_binning = os.path.join(dir_out, 'PROCESSING' ,'4_binning')


"""
CHECK INPUT FILES
"""
input_dir = config['args']['input']
# List of file paths matching the pattern
#replace R1 to 1 for SRA reads
extn=config['args']['extn']

if config['args']['sequencing'] == 'paired':
    # match only read1 files where the final pair suffix appears immediately before the extension
    pattern = os.path.join(input_dir, f'*_[R]?[1].{extn}')
    file_paths = sorted(set(glob.glob(pattern) + glob.glob(pattern)))
    r1_files = sorted(set(glob.glob(pattern)))
    
    pair_regex = re.compile(
        rf'^(?P<sample>.+?)_(?P<read_extn>R?[12])\.{re.escape(extn)}$',
        re.IGNORECASE
    )
    
    paired_reads = {}

    for fp in r1_files:
        fname = os.path.basename(fp)
        match = pair_regex.match(fname)
        if match:
            sample = match.group('sample')
            read_extn = match.group('read_extn')
            paired_reads.setdefault(sample, {})[read_extn] = fp

    for sample, reads in paired_reads.items():
        r1 = reads.get('1') or reads.get('R1')
        r2 = reads.get('2') or reads.get('R2')

        if not r1 or not r2:
            raise ValueError(f"Missing pair for sample {sample}")

        print(f"{sample}: R1={r1}, R2={r2}")
    
    #getting the number of samples 
    sample_names = list(paired_reads.keys())
    N_SAMPLES = len(sample_names)
    THRESHOLD = config['args'].get('mapping_threshold', 100)


"""ONSTART/END/ERROR
Tasks to perform at various stages the start and end of a run.
"""
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    target_log = os.path.join(dir['log'], os.path.basename(current_log))
    shutil.copy(current_log, target_log)

dir = {'log': os.path.join(dir_out, 'logs')}
onstart:
    """Cleanup old log files before starting"""
    if os.path.isdir(dir["log"]):
        oldLogs = filter(re.compile(r'.*.log').match, os.listdir(dir["log"]))
        for logfile in oldLogs:
            os.unlink(os.path.join(dir["log"], logfile))
onsuccess:
    """Print a success message"""
    sys.stderr.write('\n\n Workflow ran successfully!\n\n')

onerror:
    """Print an error message"""
    sys.stderr.write('\n\n Workflow run failed\n\n')
