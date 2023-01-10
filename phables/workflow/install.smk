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
db_files = []

db_files.append(os.path.join(DBPATH, config['phrogs_mmseqs_folder']))
db_files.append(os.path.join(DBPATH, config['smg_hmm_file']))


"""RUN SNAKEMAKE"""
rule all:
    input:
        db_files


"""RULES"""
rule phrogs_mmseqs_download:
    params:
        url=os.path.join(config['phrogs_mmseqs']),
        file=os.path.join(DBPATH, config['phrogs_mmseqs_file'])
    output:
        directory(os.path.join(DBPATH, config['phrogs_mmseqs_folder']))
    shell:
        """
            curl -Lo {params.file} {params.url}
            tar -xf {params.file} -C {DBPATH}
            rm -rf {params.file}
        """

rule smg_hmm_download:
    params:
        url=os.path.join(config['smg_hmm'])
    output:
        os.path.join(DBPATH, config['smg_hmm_file'])
    shell:
        """
            curl -Lo {output} {params.url}
        """