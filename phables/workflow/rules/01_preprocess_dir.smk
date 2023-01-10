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

INDIR = config['input']
READDIR = config['reads']
OUTDIR = config['output']
print(f"Output files will be saved to directory, {OUTDIR}\n")


############################################################################
# Checking for the assembly graph file
############################################################################

GRAPH_FILE = os.path.join(INDIR, 'processing', 'assembly', 'CONTIG_DICTIONARY', 'FLYE', 'assembly_graph.gfa')

if not os.path.exists(GRAPH_FILE):
    sys.stderr.write("ERROR: Could not find the assembly graph file from the input.\n")
    sys.exit(0)


############################################################################
# Checking through the reads folder
############################################################################

#READS_FILES = os.listdir(READDIR)




"""ONSTART/END/ERROR
Tasks to perform at various stages the start and end of a run.
"""
onsuccess:
    """Print a success message"""
    sys.stderr.write('\n\nSuccess!\n\n')

onerror:
    """Print an error message"""
    sys.stderr.write('\n\nFailed\n\n')