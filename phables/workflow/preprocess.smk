"""
Phables: Phage bubbles resolve bacteriophage genomes in viral metagenomic samples.

2023, Vijini Mallawaarachchi

This is an auxiliary Snakefile to format unitigs and preprocess.
"""


"""CONFIGURATION"""
configfile: os.path.join(workflow.basedir, 'config', 'config.yaml')
configfile: os.path.join(workflow.basedir, 'config', 'databases.yaml')


"""PREFLIGHT CHECKS
Validate your inputs, set up directories, parse your config, etc.
"""
include: "rules/01_preprocess_dir.smk"


"""TARGETS
Declare your targets, either here, or in a separate file.
"""
include: "rules/02_preprocess_targets.smk"


"""RUN SNAKEMAKE!"""
rule all:
    input:
        preprocess_files


"""RULES
Add rules files with the include directive here, or add rules AFTER rule 'all'.
"""
include: "rules/gfa2fasta.smk"