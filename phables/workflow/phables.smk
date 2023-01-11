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


# Mark target rules
target_rules = []
def targetRule(fn):
    assert fn.__name__.startswith('__')
    target_rules.append(fn.__name__[2:])
    return fn


"""RUN SNAKEMAKE!"""
@targetRule
rule all:
    input:
        allTargets


"""RULES
Add rules files with the include directive here, or add rules AFTER rule 'all'.
"""

@targetRule
rule run_phables:
    input:
        GRAPH_FILE,
        INFO_FILE,
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
    log:
        os.path.join(LOGSDIR, "phables_output.log")
    conda: 
        "./envs/phables.yaml"
    shell:
        """
            phables/workflow/scripts/phables.py -g {GRAPH_FILE} -p {INFO_FILE} -hm {SMG_FILE} -ph {PHROG_ANNOT} -c {COVERAGE_FILE} -b {BAM_PATH} -o {OUTDIR} -l {log}
        """