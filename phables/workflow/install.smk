"""
Phables: Phage bubbles resolve bacteriophage genomes in viral metagenomic samples.

2023, Vijini Mallawaarachchi

This is an auxiliary Snakefile to install databases or dependencies.
"""


"""CONFIGURATION"""
configfile: os.path.join(workflow.basedir, 'config', 'config.yaml')
configfile: os.path.join(workflow.basedir, 'config', 'databases.yaml')

include: "rules/00_database_dir.smk"


"""TARGETS"""
allDatabaseFiles = []

allDatabaseFiles.append(os.path.join(databaseDir, config['phrogs_mmseqs_file']))
allDatabaseFiles.append(os.path.join(databaseDir, config['smg_hmm_file']))


"""RUN SNAKEMAKE"""
rule all:
    input:
        allDatabaseFiles


"""RULES"""
rule phrogs_mmseqs_download:
    params:
        url=os.path.join(config['phrogs_mmseqs'])
    output:
        os.path.join(databaseDir, config['phrogs_mmseqs_file'])
    shell:
        """
            curl -Lo {output} {params.url}
        """

rule smg_hmm_download:
    params:
        url=os.path.join(config['smg_hmm'])
    output:
        os.path.join(databaseDir, config['smg_hmm_file'])
    shell:
        """
            curl -Lo {output} {params.url}
        """