"""
Use raw_coverage to map to calculate coverage of unitigs.
Use combine_cov to combine the coverage values of multiple samples into one file.
"""

rule koverage_tsv:
    """Generate TSV of samples and reads for Koverage"""
    output:
        os.path.join(OUTDIR,"samples.tsv")
    params:
        SAMPLE_READS
    run:
        from metasnek import fastq_finder
        fastq_finder.write_samples_tsv(params[0], output[0])


rule koverage:
    """Get coverage statistics with Koverage + CoverM"""
    input:
        tsv = os.path.join(OUTDIR,"samples.tsv"),
        edges = EDGES_FILE
    params:
        out_dir = OUTDIR
    output:
        expand(os.path.join(OUTDIR, "temp", "{sample}.{ext}"),
               sample=SAMPLE_NAMES,
               ext=["bam","bam.bai"]),
        os.path.join(OUTDIR, "results", "sample_bench_coverage.tsv")
    threads:
        config["resources"]["jobCPU"]
    resources:
        mem_mb = config["resources"]["jobMem"]
    conda:
        os.path.join("..", "envs", "koverage.yaml")
    shell:
        """
        koverage run bench \
            --reads {input.tsv} \
            --ref {input.edges} \
            --threads {threads} \
            --output {params.out_dir}
        """


rule run_combine_cov:
    """Sample\tContig\tCount\tRPM\tRPKM\tRPK\tTPM\tMean\tMedian\tHitrate\tVariance\n"""
    input:
        os.path.join(OUTDIR, "results", "sample_bench_coverage.tsv")
    output:
        os.path.join(OUTDIR, "coverage.tsv")
    shell:
        """awk -F '\t' '{{ sum[$2] += $4 }} END {{ for (key in sum) print key, sum[key] }}' {input} > {output}"""
