"""
Use CoverM to map to calculate coverage of unitigs.
Use combine_cov to combine the coverage values of multiple samples into one file.
"""

rule run_coverm:
    input:
        edges = EDGES_FILE,
        r1 = os.path.join(READ_DIR, PATTERN_R1),
        r2 = os.path.join(READ_DIR, PATTERN_R2)
    output:
        os.path.join(COVERM_PATH, "{sample}_rpkm.tsv")
    log:
        os.path.join(LOGSDIR, "{sample}_coverm.log")
    conda: 
        "../envs/mapping.yaml"
    shell:
        """
            TMPDIR={OUTDIR}
            coverm contig -m rpkm -1 {input.r1} -2 {input.r2} -r {input.edges} -t {threads} --output-file {output}
        """


rule run_combine_cov:
    input:
        files=expand(os.path.join(COVERM_PATH, "{sample}_rpkm.tsv"), sample=SAMPLES)
    output:
        os.path.join(OUTDIR, "coverage.tsv")
    log:
        os.path.join(LOGSDIR, "combine_cov.log")
    conda: 
        "../envs/phables.yaml"
    shell:
        """
            python phables/workflow/scripts/combine_cov.py --covpath {COVERM_PATH} --output {OUTDIR} > {log}
        """