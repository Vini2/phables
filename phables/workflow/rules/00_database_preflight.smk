"""
Add your preflight checks as pure Python code here.
e.g. Configure the run, declare directories, validate the input files etc.
This preflight check to confirm the database filepaths 
"""

"""CONFIGURATION"""
configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')
configfile: os.path.join(workflow.basedir, '../', 'config', 'databases.yaml')

"""CHECK IF CUSTOM DATABASE DIRECTORY"""
DBPATH = ""
if config['db_dir'] is None:
    DBPATH = os.path.join(workflow.basedir, 'databases')
else:
    DBPATH = config['db_dir']


"""ONSTART/END/ERROR
Tasks to perform at various stages the start and end of a run.
"""
onsuccess:
    """Print a success message"""
    sys.stderr.write('\n\nDatabases are successfully setup!\n\n')

onerror:
    """Print an error message"""
    sys.stderr.write('\n\nERROR: Databases were not setup.\n\n')