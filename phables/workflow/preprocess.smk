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
include: "rules/00_database_preflight.smk"
include: "rules/01_preprocess_preflight.smk"


"""TARGETS
Declare your targets, either here, or in a separate file.
"""
include: "rules/01_preprocess_targets.smk"


"""RUN SNAKEMAKE!"""
rule all:
    input:
        allTargets


"""RULES
Add rules files with the include directive here, or add rules AFTER rule 'all'.
"""

# Step 2: Obtain unitig sequences from assembly graph 
include: "rules/gfa2fasta.smk"

# Step 3: Map reads to unitig sequences and get BAM files
include: "rules/mapping.smk"

# Step 4: Run CoverM to get coverage of unitig sequences
include: "rules/coverm.smk"

# Step 5: Scan unitig sequences for single-copy marker genes and PHROGs
include: "rules/genes.smk"