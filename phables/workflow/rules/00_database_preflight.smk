"""
Add your preflight checks as pure Python code here.
e.g. Configure the run, declare directories, validate the input files etc.
This preflight check to confirm the database filepaths 
"""


"""CHECK IF CUSTOM DATABASE DIRECTORY"""
DBPATH = ""
try:
    if config['databases'] is None:
        DBPATH = os.path.join(workflow.basedir, '..', '..', 'databases')
    else:
        DBPATH = config['databases']
except KeyError:
    DBPATH = os.path.join(workflow.basedir,'..','..','databases')


"""ONSTART/END/ERROR
Tasks to perform at various stages the start and end of a run.
"""
onsuccess:
    """Print a success message"""
    sys.stderr.write('\n\nDatabases are successfully setup!\n\n')

onerror:
    """Print an error message"""
    sys.stderr.write('\n\nERROR: Databases were not setup.\n\n')
