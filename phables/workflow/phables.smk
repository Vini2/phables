"""
Phables: Phage bubbles resolve bacteriophage genomes in viral metagenomic samples.

2023, Vijini Mallawaarachchi

This is the main Snakefile to run phables.
"""

"""CONFIGURATION"""
configfile: os.path.join(workflow.basedir, 'config', 'config.yaml')
configfile: os.path.join(workflow.basedir, 'config', 'databases.yaml')


"""PREFLIGHT CHECKS
Validate your inputs, set up directories, parse your config, etc.
"""
include: "rules/00_database_preflight.smk"
include: "rules/02_phables_preflight.smk"


"""TARGETS
Declare your targets, either here, or in a separate file.
"""
include: "rules/02_phables_targets.smk"


"""RUN SNAKEMAKE!"""
rule all:
    input:
        allTargets


"""RULES
Add rules files with the include directive here, or add rules AFTER rule 'all'.
"""

rule run_phables:
    input:
        GRAPH_FILE,
        COVERAGE_FILE,
        BAM_PATH,
        PHROG_ANNOT,
        SMG_FILE
    output:
        genomes_fasta = os.path.join(OUTDIR, "resolved_paths.fasta"),
        genomes_folder = directory(os.path.join(OUTDIR, "resolved_phages")),
        genome_info = os.path.join(OUTDIR, "resolved_genome_info.txt"),
        unitigs = os.path.join(OUTDIR, "resolved_edges.fasta"),
        component_info = os.path.join(OUTDIR, "resolved_component_info.txt")
    params:
        graph = GRAPH_FILE,
        hmmout = SMG_FILE,
        phrogs = PHROG_ANNOT,
        coverage = COVERAGE_FILE,
        bampath = BAM_PATH,
        minlength = ML,
        mincov = MC,
        compcount = CC,
        maxpaths = MP,
        mgfrac = MGF,
        alignscore = AS,
        seqidentity = SI,
        output = OUTDIR,
        log = os.path.join(LOGSDIR, "phables_output.log")
    log:
        os.path.join(LOGSDIR, "phables_output.log")
    conda: 
        "./envs/phables.yaml"
    script:
        os.path.join('scripts', 'phables.py')