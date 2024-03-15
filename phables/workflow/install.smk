"""
Phables: from fragmented assemblies to high-quality bacteriophage genomes.

2023, Vijini Mallawaarachchi

This is an auxiliary Snakefile to install databases or dependencies.
"""


"""CONFIGURATION"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")
configfile: os.path.join(workflow.basedir, "..", "config", "databases.yaml")

include: "rules/00_database_preflight.smk"


"""TARGETS"""
db_files = []

db_files.append(os.path.join(DBPATH, config['phrogs_mmseqs_folder']))
db_files.append(os.path.join(DBPATH, config['smg_hmm_file']))
db_files.append(os.path.join(DBPATH, config['phrog_annot_file']))


"""RUN SNAKEMAKE"""
rule all:
    input:
        db_files


"""RULES"""
rule phrogs_mmseqs_download:
    params:
        url=os.path.join(config['phrogs_mmseqs']),
        file=os.path.join(DBPATH, config['phrogs_mmseqs_file']),
        db_path = DBPATH
    output:
        directory(os.path.join(DBPATH, config['phrogs_mmseqs_folder']))
    conda:
        os.path.join("envs", "curl.yaml")
    shell:
        """
            curl -Lko {params.file} {params.url}
            tar -xf {params.file} -C {params.db_path}
            rm -rf {params.file}
        """

rule smg_hmm_download:
    params:
        url=os.path.join(config['smg_hmm'])
    output:
        os.path.join(DBPATH, config['smg_hmm_file'])
    conda:
        os.path.join("envs", "curl.yaml")
    shell:
        """
            curl -Lko {output} {params.url}
        """

rule phrog_annot_download:
    params:
        url=os.path.join(config['phrog_annot'])
    output:
        os.path.join(DBPATH, config['phrog_annot_file'])
    conda:
        os.path.join("envs", "curl.yaml")
    shell:
        """
            curl -Lko {output} {params.url}
        """