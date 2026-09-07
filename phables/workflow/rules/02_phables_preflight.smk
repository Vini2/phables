"""
Add your preflight checks as pure Python code here.
e.g. Configure the run, declare directories, validate the input files etc.
This preflight check to confirm the database filepaths 
"""

from metasnek import fastq_finder

"""
Setting the directory variables
"""

# THREADS = config['threads']
INPUT = config['input']
OUTDIR = config['output']
print(f"Output files will be saved to directory, {OUTDIR}\n")


############################################################################
# Checking through the reads folder
############################################################################

SAMPLE_READS = fastq_finder.parse_samples_to_dictionary(config['reads'])
SAMPLE_NAMES = list(SAMPLE_READS.keys())


############################################################################
# Get Phables parameters
############################################################################
ML = config['minlength']
MC = config['mincov']
CC = config['compcount']
MP = config['maxpaths']
# .get(): both are optional and absent from older config files.
MFD_TIME_LIMIT = config.get('mfd_time_limit')
MFD_DUMP_SLOW = config.get('mfd_dump_slow')
MGF = config['mgfrac']
GC = config['genecaller']
PD = config['phagedetection']
GPU_BACKEND = config['gpu_backend']
FOLDSEEK_GPU = config['foldseek_gpu']
MFD_WORKERS = config['mfd_workers']

# Per-rule CPU/memory request. --job-cpu / --job-mem override resources.jobCPU /
# resources.jobMem when given; otherwise the config value stands.
#
# These exist because --threads does NOT do what it looks like. --threads sets
# Snakemake's total core budget (--cores); each rule separately asks for
# resources.jobCPU, which defaults to 8. Snakemake then gives the rule
# min(jobCPU, cores). So `--threads 64` on a 64-core node runs every rule 8-wide
# and leaves 56 cores idle -- measured on a real run: coverm_map and
# scan_hallmark both took their full time on 8 threads of a 64-core allocation.
#
# `or` rather than a dict .get default: click passes None when the flag is
# absent, and None is exactly the value that must fall through to the config.
# The click defaults are None for the same reason -- snaketool merges CLI
# defaults OVER config files, so any non-None default would make
# resources.jobCPU permanently unsettable from a config file.
JOB_CPU = config.get('job_cpu') or config['resources']['jobCPU']
JOB_MEM = config.get('job_mem') or config['resources']['jobMem']
EV = config['evalue']
SI = config['seqidentity']
CT = config['covtol']
AL = config['alpha']
LR = config['longreads']
PR = config['prefix']


"""DIRECTORIES/FILES etc.
Declare some directories for pipeline intermediates and outputs.
"""
LOGSDIR = os.path.join(OUTDIR, 'logs')


"""ONSTART/END/ERROR
Tasks to perform at various stages the start and end of a run.
"""
onstart:
    """Cleanup old log files before starting"""
    if os.path.isdir(LOGSDIR):
        oldLogs = filter(re.compile(r'.*.log').match, os.listdir(LOGSDIR))
        for logfile in oldLogs:
            os.unlink(os.path.join(LOGSDIR, logfile))


onsuccess:
    """Print a success message"""
    sys.stderr.write('\n\nPhables ran successfully!\n\n')


onerror:
    """Print an error message"""
    sys.stderr.write('\n\nPhables run failed\n\n')
