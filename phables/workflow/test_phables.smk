"""
Phables: Phage bubbles resolve bacteriophage genomes in viral metagenomic samples.

2023, Vijini Mallawaarachchi

This is an auxiliary Snakefile to test phables.
"""

"""CONFIGURATION"""
configfile: os.path.join(workflow.basedir, 'config', 'config.yaml')
configfile: os.path.join(workflow.basedir, 'config', 'databases.yaml')


"""PREFLIGHT CHECKS
Validate your inputs, set up directories, parse your config, etc.
"""
include: "rules/00_database_preflight.smk"
include: "rules/03_test_preflight.smk"


"""TARGETS
Declare your targets, either here, or in a separate file.
"""
include: "rules/03_test_targets.smk"


"""RUN SNAKEMAKE!"""
rule all:
    input:
        allTargets


"""RULES
Add rules files with the include directive here, or add rules AFTER rule 'all'.
"""

rule test_phables:
    input:
        g = os.path.join(TESTDIR, "assembly_graph.gfa"),
        c = os.path.join(TESTDIR, "edge_coverages.tsv"),
        b = TESTDIR,
        ph = os.path.join(TESTDIR, "phrogs_annotations.tsv"),
        hm = os.path.join(TESTDIR, "edges.fasta.hmmout")
    output:
        genomes_fasta = temp(os.path.join(TESTDIR, "resolved_paths.fasta")),
        genomes_folder = temp(directory(os.path.join(TESTDIR, "resolved_phages"))),
        genome_info = temp(os.path.join(TESTDIR, "resolved_genome_info.txt")),
        phage_edges = temp(os.path.join(TESTDIR, "phage_like_edges.fasta")),
        unitigs = temp(os.path.join(TESTDIR, "resolved_edges.fasta")),
        component_info = temp(os.path.join(TESTDIR, "resolved_component_info.txt")),
        mfd = temp(os.path.join(TESTDIR, "results_MFD.txt")),
        mfd_details = temp(os.path.join(TESTDIR, "results_MFD_details.txt")),
        log = temp(os.path.join(TESTDIR, "phables_output.log"))
    params:
        graph = os.path.join(TESTDIR, "assembly_graph.gfa"),
        hmmout = os.path.join(TESTDIR, "edges.fasta.hmmout"),
        phrogs = os.path.join(TESTDIR, "phrogs_annotations.tsv"),
        coverage = os.path.join(TESTDIR, "edge_coverages.tsv"),
        bampath = TESTDIR,
        minlength = ML,
        mincov = MC,
        compcount = CC,
        maxpaths = MP,
        mgfrac = MGF,
        alignscore = AS,
        seqidentity = SI,
        output = TESTDIR,
        log = temp(os.path.join(TESTDIR, "phables_output.log"))
    log:
        temp(os.path.join(TESTDIR, "phables_output.log"))
    conda: 
        "./envs/phables.yaml"
    script:
        os.path.join('scripts', 'phables.py')