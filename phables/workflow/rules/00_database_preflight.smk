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


"""HALLMARK DB DEFAULTS
--hallmark-db/--hallmark-categories default to None in config.yaml, so a plain
run doesn't need them unless --phagedetection prostt5-foldseek is actually
used (see genes.smk's own note on why this has to stay gated rather than an
always-present path/input declaration). Left as None, they'd previously have
to be passed explicitly every time. Now that `phables install` fetches a
pre-built hallmark_db (install.smk's hallmark_db_download rule) under DBPATH,
resolve to that location automatically when unset -- same "just works once
installed" behaviour --databases's own marker.hmm/phrogs_mmseqs_db already
have -- while an explicit --hallmark-db/--hallmark-categories (e.g. pointing
at your own rebuild) still overrides this, since we only fill in a None.
Harmless to run in install.smk/test_phables.smk too (this file is included by
all three) -- neither of those reads hallmark_db/hallmark_categories, so
setting the config value there does nothing.
"""
if config.get('hallmark_db') is None:
    config['hallmark_db'] = os.path.join(DBPATH, config['hallmark_db_folder'], 'hallmark_db')
if config.get('hallmark_categories') is None:
    config['hallmark_categories'] = os.path.join(DBPATH, config['hallmark_db_folder'], 'hallmark_categories.tsv')


"""ONSTART/END/ERROR
Tasks to perform at various stages the start and end of a run.
"""
onsuccess:
    """Print a success message"""
    sys.stderr.write('\n\nDatabases are successfully setup!\n\n')

onerror:
    """Print an error message"""
    sys.stderr.write('\n\nERROR: Databases were not setup.\n\n')
