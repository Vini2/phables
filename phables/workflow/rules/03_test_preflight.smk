"""
Add your preflight checks as pure Python code here.
e.g. Configure the run, declare directories, validate the input files etc.
This preflight check to confirm the database filepaths 
"""


"""
Setting the directory variables
"""

TESTDIR = config['dir']


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


"""ONSTART/END/ERROR
Tasks to perform at various stages the start and end of a run.
"""
onsuccess:
    """Print a success message"""
    sys.stderr.write('\n\nPhables test run was successful!\n\n')

onerror:
    """Print an error message"""
    sys.stderr.write('\n\nPhables test run failed! Please check.\n\n')