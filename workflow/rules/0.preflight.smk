#Snakemake rule to amke sure the input files are correct and all the params are provided to start the workflow

#To start with the input files have to be in a folders
#can be compressed or decompressed
#For now there can be no sub folders, the files have to be saved in the same folder

import os 
import glob
import yaml
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
dir_assembly = os.path.join(dir_out, 'PROCESSING' ,'3_coassembly')

"""
CHECK INPUT FILES
"""
input_dir = config['args']['input']
# List of file paths matching the pattern
#replace R1 to 1 for SRA reads
extn=config['args']['extn']

if config['args']['sequencing'] == 'paired':
    patten= os.path.join(input_dir, f'*_1*.{extn}')
    file_paths = glob.glob(patten)
    samples_names = [os.path.splitext(os.path.basename(file_path))[0].rsplit('_1', 1)[0] for file_path in file_paths]
elif config['args']['sequencing'] == 'longread':
    pattern = os.path.join(input_dir, f"*.{extn}")
    file_paths = glob.glob(pattern)
    samples_names = (os.path.splitext(os.path.basename(file_path)) if '.' in os.path.basename(file_path) else (os.path.basename(file_path), ''))

print(f"Samples are {samples_names}")

FQEXTN = extn[0]
PATTERN_R1 = '{sample}_1.' + extn
PATTERN_R2 = '{sample}_2.' + extn
PATTERN_LONG= '{sample}' + extn

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
