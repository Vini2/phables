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
