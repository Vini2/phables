"""
Add your preflight checks as pure Python code here.
e.g. Configure the run, declare directories, validate the input files etc.
This preflight check to confirm the database filepaths 
"""


"""
Setting the directory variables
"""

THREADS = config['threads']
INPUT = config['input']
OUTDIR = config['output']
print(f"Output files will be saved to directory, {OUTDIR}\n")


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


############################################################################
# Checking for required input files for phables
############################################################################

# Check for assembly_graph.gfa
GRAPH_FILE = INPUT
if not os.path.exists(GRAPH_FILE):
    sys.stderr.write("ERROR: Could not find the assembly_graph.gfa file from the input.\n")
    sys.exit(0)

# Check for unitigs file
EDGES_FILE = os.path.join(OUTDIR, 'edges.fasta')
if not os.path.exists(EDGES_FILE):
    sys.stderr.write("ERROR: Could not find the edges.fasta file. Please make sure to run the preprocessing steps.\n")
    sys.exit(0)

# Check for coverage file
COVERAGE_FILE = os.path.join(OUTDIR, 'coverage.tsv')
if not os.path.exists(EDGES_FILE):
    sys.stderr.write("ERROR: Could not find the coverage.tsv file. Please make sure to run the preprocessing steps.\n")
    sys.exit(0)

# Check for BAM folder
BAM_PATH = os.path.join(OUTDIR, 'bam_files')
if not os.path.exists(BAM_PATH):
    sys.stderr.write("ERROR: Could not find the path to BAM files. Please make sure to run the preprocessing steps.\n")
    sys.exit(0)

# Check if there are BAM files in BAM folder
SAMPLES = glob_wildcards(os.path.join(BAM_PATH, '{sample}.bam'))
if len(SAMPLES) == 0:
    sys.stderr.write("ERROR: Could not find any BAM files in {BAM_PATH}. Please make sure to run the preprocessing steps.\n")
    sys.exit(0)

# Check if there are index files of BAM files in BAM folder
SAMPLES = glob_wildcards(os.path.join(BAM_PATH, '{sample}.bam.bai'))
if len(SAMPLES) == 0:
    sys.stderr.write("ERROR: Could not find any index files for the BAM files in {BAM_PATH}. Please make sure to run the preprocessing steps.\n")
    sys.exit(0)

# Check for PHROG annotations
PHROGS_PATH = os.path.join(OUTDIR, 'phrogs/')
PHROG_ANNOT = os.path.join(PHROGS_PATH, 'phrogs_annotations.tsv')
if not os.path.exists(PHROG_ANNOT):
    sys.stderr.write("ERROR: Could not find the PHROG annotations. Please make sure to run the preprocessing steps.\n")
    sys.exit(0)

# Check for single-copy marker gene results file
SMG_FILE = os.path.join(OUTDIR, 'edges.fasta.hmmout')
if not os.path.exists(SMG_FILE):
    sys.stderr.write("ERROR: Could not find the bacterial single-copy marker gene annotations. Please make sure to run the preprocessing steps.\n")
    sys.exit(0)


############################################################################
# Get Phables parameters
############################################################################
ML = config['minlength']
MC = config['mincov']
CC = config['compcount']
MP = config['maxpaths']
MGF = config['mgfrac']
AS = config['alignscore']
SI = config['seqidentity']


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
