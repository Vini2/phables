"""
Add your preflight checks as pure Python code here.
e.g. Configure the run, declare directories, validate the input files etc.
This preflight check to confirm the database filepaths 
"""

"""CONFIGURATION"""
configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')
configfile: os.path.join(workflow.basedir, '../', 'config', 'databases.yaml')

"""
Setting the directory variables
"""

THREADS = config['threads']
INPUT = config['input']
OUTDIR = config['output']
print(f"Output files will be saved to directory, {OUTDIR}\n")


############################################################################
# Checking for the assembly graph files
############################################################################

GRAPH_FILE = INPUT
if not os.path.exists(GRAPH_FILE):
    sys.stderr.write("ERROR: Could not find the .GFA file.\n")
    sys.exit(0)


############################################################################
# Checking through the reads folder
############################################################################

READ_DIR = config['reads']
SAMPLES,EXTENSIONS, = glob_wildcards(os.path.join(READ_DIR, '{sample}_R1{extn}'))

# Check if there are read files
if len(SAMPLES) == 0:
    sys.stderr.write("ERROR: Could not find any FASTQ files in {READ_DIR}. Please check the reads path.\n")
    sys.exit(0)

if len(set(EXTENSIONS)) != 1:
    sys.stderr.write("ERROR: You have more than one type of file extension. Please make sure that you have the same file extension.\n")
    sys.exit(0)

FQEXTN = EXTENSIONS[0]
PATTERN_R1 = '{sample}_R1' + FQEXTN
PATTERN_R2 = '{sample}_R2' + FQEXTN


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


"""ONSTART/END/ERROR
Tasks to perform at various stages the start and end of a run.
"""
onsuccess:
    """Print a success message"""
    sys.stderr.write('\n\nDirectory preprocessing success!\n\n')

onerror:
    """Print an error message"""
    sys.stderr.write('\n\nDirectory preprocessing failed\n\n')