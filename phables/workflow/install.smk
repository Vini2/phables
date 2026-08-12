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
db_files.append(os.path.join(DBPATH, config['hallmark_db_folder']))


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

rule hallmark_db_download:
    """See docs/hallmark_db.md for what this is and how it's built
    (workflow/scripts/build_hallmark_db.py) -- this rule only fetches the
    pre-built tarball, it doesn't build anything from phold's own DB.

    Expects the tarball's members at its OWN root (i.e. built with
    `cd hallmark_db/ && tar -czhf hallmark_db.tar.gz .` -- note -h/
    --dereference, needed because createsubdb leaves .lookup/.source as
    symlinks back to the original phold DB's own files; skipping -h ships
    broken links to everyone except the machine it was built on), NOT nested
    inside its own hallmark_db/ directory -- unlike phrogs_mmseqs_download
    above, which extracts straight to DBPATH because that tarball already
    contains its own top-level folder. Extracting this one the same way would
    double-nest it (DBPATH/hallmark_db/hallmark_db/hallmark_db...), so this
    creates {output} itself first and extracts into that instead.
    """
    params:
        url=os.path.join(config['hallmark_db_url']),
        file=os.path.join(DBPATH, config['hallmark_db_file']),
    output:
        directory(os.path.join(DBPATH, config['hallmark_db_folder']))
    conda:
        os.path.join("envs", "curl.yaml")
    shell:
        """
            curl -Lko {params.file} {params.url}
            mkdir -p {output}
            tar -xf {params.file} -C {output}
            rm -rf {params.file}
        """