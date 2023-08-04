rule koverage_genomes:
    """Get coverage statistics with Koverage"""
    input:
        tsv = os.path.join(OUTDIR,"phables.samples.tsv"),
        genomes = RESOLVED_GENOMES
    params:
        out_dir = OUTDIR,
        profile = lambda wildcards: "--profile " + config["profile"] if config["profile"] else "",
    output:
        os.path.join(OUTDIR, "results", "sample_coverage.tsv")
    threads:
        config["resources"]["jobCPU"]
    resources:
        mem_mb = config["resources"]["jobMem"],
        mem = str(config["resources"]["jobMem"]) + "MB"
    conda:
        os.path.join("..", "envs", "koverage.yaml")
    shell:
        """
        koverage run \
            --reads {input.tsv} \
            --ref {input.genomes} \
            --threads {threads} \
            --output {params.out_dir} \
            {params.profile}
        """
