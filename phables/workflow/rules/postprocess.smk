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


"""
Per-genome coverage of the FINAL resolved genomes (+ unresolved phage-like
edges), for phables' human-facing report tables (sample_genome_read_counts.tsv,
sample_genome_rpkm.tsv, sample_genome_mean_coverage.tsv). This does NOT feed the
MFD flow decomposition -- that uses preprocess/coverage.tsv, built separately in
coverage.smk.

The three rules below replace the former single `koverage run` call, completing
PLAN.md §4.8 ("Drop the koverage wrapper -- call CoverM directly"). coverage.smk
converted the preprocess step first; this rule was the LAST remaining thing in
the whole workflow that actually depended on koverage (the other two
koverage-prefixed rules were misnomers -- `samples_tsv` in coverage.smk uses
metasnek, and `format_genome_coverage` below uses phables.yaml; neither ever
invoked koverage). With this converted, koverage.yaml is deleted and §4.8's
headline benefit -- "removes a nested Snakemake dependency from inside what is
already a Snakemake workflow" -- is finally actually realised, which the
preprocess-only change had NOT achieved on its own.

Structurally identical to coverage.smk's coverm_map / coverm_bam2counts /
coverm_combine, just against genomes_and_unresolved_edges.fasta instead of
edges.fasta, writing under postprocess/ instead of preprocess/. Kept as a
parallel copy rather than shared/parameterised because Snakemake rules are
declarative -- the two differ only in input/output paths.

IMPORTANT -- the estimator changed here, unlike in coverage.smk. The old
`koverage run` (no `coverm` subcommand) used Koverage's *native* "map" mode,
whose columns are Sample/Contig/Count/RPM/RPKM/RPK/TPM/Mean/Median/Hitrate/
Variance. `coverm contig` emits Sample/Contig/Count/RPKM/TPM/Mean/
Covered_fraction/Variance instead. The three values phables actually consumes
(Count, RPKM, Mean) all still exist, and format_koverage_results.py's column
indices are updated to match -- but native-koverage and CoverM compute them via
different code paths (Koverage's own minimap2+python vs CoverM over a sorted
BAM), so the reported numbers may differ slightly. That is a REPORT-ONLY
difference (these tables are not read by anything downstream), but it has not
been numerically regression-tested against the old output on real data -- worth
doing per PLAN.md §4.8's own advice before relying on the absolute values.
"""

rule coverm_map_genomes:
    """Map each sample's reads to the resolved genomes -> sorted,
    unmapped-filtered, indexed BAM. Mirrors coverage.smk's coverm_map."""
    input:
        ref = os.path.join(OUTDIR, "postprocess", "genomes_and_unresolved_edges.fasta"),
        r1 = lambda wildcards: SAMPLE_READS[wildcards.sample]["R1"],
    params:
        r2 = lambda wildcards: SAMPLE_READS[wildcards.sample]["R2"] or "",
        preset = "map-ont" if LR else "sr",
    output:
        # temp(): unlike the preprocess BAMs (which 02_phables_targets.smk
        # declares as real targets because coverage_utils.py's pysam mate-pair
        # /junction logic reads them back), nothing consumes these after
        # coverm_bam2counts_genomes below. They are also NEW artifacts -- the
        # koverage native "map" mode this replaced piped minimap2 straight into
        # a python counter and never wrote a postprocess BAM at all -- so
        # keeping them would be a pure disk regression, which matters at
        # 10k-sample scale. Use Snakemake's --notemp to retain them for
        # debugging.
        bam = temp(os.path.join(OUTDIR, "postprocess", "temp", "{sample}.bam")),
        bai = temp(os.path.join(OUTDIR, "postprocess", "temp", "{sample}.bam.bai")),
    threads:
        config["resources"]["jobCPU"]
    resources:
        mem_mb = config["resources"]["jobMem"]
    conda:
        None if CONTAINER_IMAGE else os.path.join("..", "envs", "coverm.yaml")
    container:
        CONTAINER_IMAGE
    log:
        os.path.join(LOGSDIR, "coverm_map_genomes.{sample}.log")
    shell:
        """
        {{ minimap2 -t {threads} -ax {params.preset} --secondary=no {input.ref} {input.r1} {params.r2} \
            | samtools sort -T {wildcards.sample}_genomes -@ {threads} - \
            | samtools view -b -F 4 > {output.bam} ; \
        samtools index {output.bam} ; }} 2> {log}
        """


rule coverm_bam2counts_genomes:
    """Per-sample coverage stats over the resolved genomes."""
    input:
        os.path.join(OUTDIR, "postprocess", "temp", "{sample}.bam")
    output:
        temp(os.path.join(OUTDIR, "postprocess", "temp", "{sample}.cov"))
    conda:
        None if CONTAINER_IMAGE else os.path.join("..", "envs", "coverm.yaml")
    container:
        CONTAINER_IMAGE
    log:
        os.path.join(LOGSDIR, "coverm_bam2counts_genomes.{sample}.log")
    shell:
        """
        coverm contig -b {input} \
            -m count -m rpkm -m tpm -m mean -m covered_fraction -m variance \
            > {output} 2> {log}
        """


rule coverm_combine_genomes:
    """Reshape per-sample coverm TSVs into one long-format
    Sample/Contig/<method...> table. Mirrors coverage.smk's coverm_combine.
    Output keeps its original sample_coverage.tsv path/name so the report step
    below is unchanged, but its COLUMNS are now coverm-mode, not koverage
    native-map-mode -- see this file's header."""
    input:
        expand(os.path.join(OUTDIR, "postprocess", "temp", "{sample}.cov"), sample=SAMPLE_NAMES)
    output:
        os.path.join(OUTDIR, "postprocess", "results", "sample_coverage.tsv")
    run:
        with open(input[0]) as f:
            header = f.readline().rstrip("\n").split("\t")
        header = [" ".join(col.split()[1:]) for col in header]
        header[0] = "Contig"
        with open(output[0], "w") as out:
            out.write("Sample\t" + "\t".join(header) + "\n")
            for sample, cov_file in zip(SAMPLE_NAMES, input):
                with open(cov_file) as f:
                    next(f)  # skip this sample's own header line
                    for line in f:
                        out.write(f"{sample}\t{line}")


rule format_genome_coverage:
    """Format the per-genome coverage table into phables' report TSVs."""
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
        None if CONTAINER_IMAGE else os.path.join("..", "envs", "phables.yaml")
    container:
        CONTAINER_IMAGE
    script:
        os.path.join("..", "scripts", "format_koverage_results.py")