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
MGF = config['mgfrac']
GC = config['genecaller']
PD = config['phagedetection']
GPU_BACKEND = config['gpu_backend']
FOLDSEEK_GPU = config['foldseek_gpu']
# CONTAINER_IMAGE: workflow-wide container (e.g. a phables release image with
# every per-rule tool baked in), replacing --use-conda's per-rule env creation
# entirely when set -- every rule below does
# `conda: None if CONTAINER_IMAGE else os.path.join(...)` +
# `container: CONTAINER_IMAGE`, so a real value here disables conda for ALL of
# them at once (not just predict_3di). Run with --use-singularity, no
# --use-conda, once this is set.
#
# PROSTT5_CONTAINER falls back to CONTAINER_IMAGE when not set explicitly --
# --prostt5-container remains available to point JUST predict_3di at a
# different image (e.g. phold's own) than the rest of the workflow, but
# --container alone is now sufficient to route everything, predict_3di
# included, through one image.
PROSTT5_CONTAINER = config['prostt5_container'] or config['container']
CONTAINER_IMAGE = config['container']
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
