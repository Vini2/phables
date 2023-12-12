rule combine_genomes_and_unresolved_edges:
    """Combine resolved genomes and unresolved edges"""
    input:
        genomes = RESOLVED_GENOMES,
        unresolved_edges = os.path.join(OUTDIR, "phables", "unresolved_phage_like_edges.fasta")
    output:
        os.path.join(OUTDIR, "postprocess", "genomes_and_unresolved_edges.fasta")
    shell:
        """
        cat {input.genomes} {input.unresolved_edges} > {output}
        """


rule koverage_genomes:
    """Get coverage statistics with Koverage"""
    input:
        tsv = os.path.join(OUTDIR, "preprocess", "phables.samples.tsv"),
        sequences = os.path.join(OUTDIR, "postprocess", "genomes_and_unresolved_edges.fasta")
    params:
        out_dir = os.path.join(OUTDIR, "postprocess"),
        profile = lambda wildcards: "--profile " + config["profile"] if config["profile"] else "",
    output:
        os.path.join(OUTDIR, "postprocess", "results", "sample_coverage.tsv")
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
            --no-report \
            --reads {input.tsv} \
            --ref {input.sequences} \
            --threads {threads} \
            --output {params.out_dir} \
            {params.profile}
        """


rule koverage_postprocess:
    """Format TSV of samples and reads from Koverage"""
    input:
        koverage_tsv = os.path.join(OUTDIR, "postprocess", "results", "sample_coverage.tsv"),
        samples_file = os.path.join(OUTDIR, "preprocess", "phables.samples.tsv"),
        seq_file = os.path.join(OUTDIR, "postprocess", "genomes_and_unresolved_edges.fasta")
    output:
        os.path.join(OUTDIR, "postprocess", "sample_genome_read_counts.tsv")
    params:
        koverage_tsv = os.path.join(OUTDIR, "postprocess", "results", "sample_coverage.tsv"),
        samples_file = os.path.join(OUTDIR, "preprocess", "phables.samples.tsv"),
        seq_file = os.path.join(OUTDIR, "postprocess", "genomes_and_unresolved_edges.fasta"),
        info_file = os.path.join(OUTDIR, "postprocess", "genomes_and_unresolved_edges_info.tsv"),
        output_path = os.path.join(OUTDIR, "postprocess"),
        log = os.path.join(LOGSDIR, "format_koverage_results_output.log")
    log:
        os.path.join(LOGSDIR, "format_koverage_results_output.log")
    conda:
        os.path.join("..", "envs", "phables.yaml")
    script:
        os.path.join("..", "scripts", "format_koverage_results.py")