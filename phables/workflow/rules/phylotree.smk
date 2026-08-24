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
        os.path.join("..", "envs", "phylotree.yaml")
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
        os.path.join("..", "envs", "phylotree.yaml")
    script:
        os.path.join("..", "scripts", "phylotree.py")