rule build_msa:
    input:
        RESOLVED_GENOMES,
    output:
        ALIGNED_GENOMES,
    params:
        nthreads = config["resources"]["jobCPU"],
        log = os.path.join(LOGSDIR, "mafft_output.log")
    threads:
        config["resources"]["jobCPU"]
    log:
        os.path.join(LOGSDIR, "mafft_output.log")
    conda:
        None if CONTAINER_IMAGE else os.path.join("..", "envs", "phylotree.yaml")
    container:
        CONTAINER_IMAGE
    shell:
        """
        mafft --auto --thread {threads} {input} > {output}
        """


rule build_tree:
    input:
        ALIGNED_GENOMES,
    output:
        TREE_FILE,
    params:
        aligned = ALIGNED_GENOMES,
        output = TREE_FILE,
        seed = 1,
        nthreads = config["resources"]["jobCPU"],
        log = os.path.join(LOGSDIR, "piqtree_output.log")
    threads:
        config["resources"]["jobCPU"]
    log:
        os.path.join(LOGSDIR, "piqtree_output.log")
    conda:
        None if CONTAINER_IMAGE else os.path.join("..", "envs", "phylotree.yaml")
    container:
        CONTAINER_IMAGE
    script:
        os.path.join("..", "scripts", "phylotree.py")