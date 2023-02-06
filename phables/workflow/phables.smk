"""
Phables: Phage bubbles resolve bacteriophage genomes in viral metagenomic samples.

2023, Vijini Mallawaarachchi

This is the main Snakefile to run phables.
"""

"""CONFIGURATION"""
configfile: os.path.join(workflow.basedir, 'config', 'config.yaml')
configfile: os.path.join(workflow.basedir, 'config', 'databases.yaml')


"""PREFLIGHT CHECKS"""
include: os.path.join("rules", "00_database_preflight.smk")
include: os.path.join("rules", "02_phables_preflight.smk")


"""TARGETS"""
include: os.path.join("rules", "02_phables_targets.smk")


"""Target rules"""
target_rules = []

def targetRule(fn):
    assert fn.__name__.startswith("__")
    target_rules.append(fn.__name__[2:])
    return fn

localrules: all, preprocess, phables, print_stages


"""Run stages"""
@targetRule
rule all:
    input:
        preprocessTargets,
        phablesTargets


@targetRule
rule preprocess:
    input:
        preprocessTargets


@targetRule
rule phables:
    input:
        phablesTargets


@targetRule
rule print_stages:
    run:
        print("\nIndividual Phables stages to run: \n", file=sys.stderr)
        print("* " + "\n* ".join(target_rules) + "\n\n", file=sys.stderr)


"""RULES"""
# Step 2: Obtain unitig sequences from assembly graph
include: os.path.join("rules", "gfa2fasta.smk")


# Step 3: Map reads to unitig sequences and get BAM files
include: os.path.join("rules", "mapping.smk")


# Step 4: Run CoverM to get coverage of unitig sequences
include: os.path.join("rules", "coverm.smk")


# Step 5: Scan unitig sequences for single-copy marker genes and PHROGs
include: os.path.join("rules", "genes.smk")


# Step 6: Run Phables
include: os.path.join("rules", "phables.smk")
