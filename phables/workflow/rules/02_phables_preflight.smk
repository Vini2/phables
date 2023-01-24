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

INPUT = config['input']
OUTDIR = config['output']
print(f"Output files will be saved to directory, {OUTDIR}\n")


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
onsuccess:
    """Print a success message"""
    sys.stderr.write('\n\nPhables ran successfully!\n\n')

onerror:
    """Print an error message"""
    sys.stderr.write('\n\nPhables run failed\n\n')