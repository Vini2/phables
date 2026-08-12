"""
Use raw_coverage to map to calculate coverage of unitigs.
Use combine_cov to combine the coverage values of multiple samples into one file.

The three rules below (coverm_map, coverm_bam2counts, coverm_combine) replace the
former single `koverage run coverm` call (PLAN.md §4.8, Vijini's recommendation:
"Drop the koverage wrapper — call CoverM directly"). They reproduce Koverage's own
real "coverm" mode chain exactly -- confirmed against its actual source
(github.com/beardymcjohnface/Koverage, workflow/rules/coverm.smk and its shipped
default config.yaml), not reimplemented from guesswork:
  - same minimap2 mapping command (-ax sr --secondary=no) and samtools
    sort/view -F 4/index chain (coverm_map_pe there),
  - same `coverm contig` --methods list and order: count, rpkm, tpm, mean,
    covered_fraction, variance (Koverage's own shipped default, which phables
    never overrode -- confirmed no `--profile`/coverm-params override exists
    anywhere in this repo),
  - same per-sample -> long-format Sample/Contig/<method...> TSV reshape
    (coverm_combine there).
Byte-for-byte the same `sample_coverm_coverage.tsv` shape as before, so
run_combine_cov below and coverage_utils.py's BAM globbing (both downstream
consumers) need no changes.

This file converted the preprocess (MFD-critical) step. postprocess.smk's
per-genome report coverage -- which used Koverage's separate *native*
(non-CoverM) engine -- was converted separately, in a follow-up; koverage is now
gone from the workflow entirely and koverage.yaml has been deleted.

Note: the phables-side docstring on run_combine_cov used to say "Covered_bases"
for column 7 -- that's actually Covered_fraction (Koverage's real 6th --methods
value, confirmed above). Harmless either way: only column 6 (Mean) is consumed
by the awk below, per notes/phables_audit.md §3.2.

Note: unlike Koverage's own rule (which always maps with -ax sr regardless of
read type -- confirmed against its source, a pre-existing gap this file no
longer inherits), coverm_map below branches on `LR` (config["longreads"],
set in 02_phables_preflight.smk and already used the same way -- a plain
module-level global, not a wildcard -- by phables.smk's own `longreads=LR`
param) to pick minimap2's -ax sr preset for short reads (unchanged) or
-ax map-ont for long reads. --longreads is a bare boolean CLI flag with no
ONT/PacBio distinction anywhere in phables (confirmed: no such option in
__main__.py), so a single long-read preset is what there's room to pick.

Fix History (dsmk-2026-08-10)
-----------------------------
Confirmed against ../metagenomic_phage_discovery/pilot10_longreads that every
sample actually run through this fork's --longreads path so far is Oxford
Nanopore (manifest.tsv: instrument_platform OXFORD_NANOPORE for all 10 runs),
assembled with myloasm's default (non-hifi) mode, which itself targets
Nanopore R10.4-class reads (see that project's run_pilot_assembly.sh) -- i.e.
map-ont, not map-pb/map-hifi, is the correct default for this path today, not
just a guess. run_pilot_phables.sh had flagged this exact gap as unverified
("whether koverage correctly selects minimap2 for single-end long reads")
before this fix. --secondary=no is kept for both presets (it suppresses
secondary alignments regardless of preset and coverage counting wants that
either way, not something specific to the sr preset).
"""

rule samples_tsv:
    """Generate TSV of samples and reads. Despite its former name
    (koverage_tsv) this rule never invoked koverage -- it's a plain metasnek
    fastq_finder call. Still produced because format_genome_coverage
    (postprocess.smk) reads it to get the sample-name column order for its
    report tables; not needed by coverm_map below, which uses SAMPLE_READS
    directly."""
    output:
        os.path.join(OUTDIR, "preprocess", "phables.samples.tsv")
    params:
        SAMPLE_READS
    run:
        from metasnek import fastq_finder
        fastq_finder.write_samples_tsv(params[0], output[0])


rule coverm_map:
    """Map each sample's reads to the unitig edges with minimap2 -> sorted,
    unmapped-filtered, indexed BAM. Matches Koverage's own coverm_map_pe rule
    for short reads; branches to a long-read minimap2 preset when LR is set
    (see module docstring's Fix History)."""
    input:
        ref = EDGES_FILE,
        r1 = lambda wildcards: SAMPLE_READS[wildcards.sample]["R1"],
    params:
        # Koverage's own rule passes "" for single-end samples (R2 None) --
        # minimap2 then maps r1 alone. Preserved as-is.
        r2 = lambda wildcards: SAMPLE_READS[wildcards.sample]["R2"] or "",
        # LR is a plain module-level global (config["longreads"], set once
        # for the whole run in 02_phables_preflight.smk) -- not a per-sample
        # wildcard -- so it's fine to resolve this at parse time rather than
        # via a lambda, same as phables.smk's own longreads=LR param.
        preset = "map-ont" if LR else "sr",
    output:
        bam = os.path.join(OUTDIR, "preprocess", "temp", "{sample}.bam"),
        bai = os.path.join(OUTDIR, "preprocess", "temp", "{sample}.bam.bai"),
    threads:
        config["resources"]["jobCPU"]
    resources:
        mem_mb = config["resources"]["jobMem"]
    conda:
        None if CONTAINER_IMAGE else os.path.join("..", "envs", "coverm.yaml")
    container:
        CONTAINER_IMAGE
    log:
        os.path.join(LOGSDIR, "coverm_map.{sample}.log")
    shell:
        """
        {{ minimap2 -t {threads} -ax {params.preset} --secondary=no {input.ref} {input.r1} {params.r2} \
            | samtools sort -T {wildcards.sample} -@ {threads} - \
            | samtools view -b -F 4 > {output.bam} ; \
        samtools index {output.bam} ; }} 2> {log}
        """


rule coverm_bam2counts:
    """Per-sample coverage stats. Koverage's own rule doesn't pass a thread
    count to coverm either (only Snakemake's own scheduling uses `threads:`
    here) -- preserved as-is rather than adding an unverified -t flag."""
    input:
        os.path.join(OUTDIR, "preprocess", "temp", "{sample}.bam")
    output:
        os.path.join(OUTDIR, "preprocess", "temp", "{sample}.cov")
    conda:
        None if CONTAINER_IMAGE else os.path.join("..", "envs", "coverm.yaml")
    container:
        CONTAINER_IMAGE
    log:
        os.path.join(LOGSDIR, "coverm_bam2counts.{sample}.log")
    shell:
        """
        coverm contig -b {input} \
            -m count -m rpkm -m tpm -m mean -m covered_fraction -m variance \
            > {output} 2> {log}
        """


rule coverm_combine:
    """Reshape per-sample coverm TSVs (each header column is
    "<bam filename> <method>") into one long-format Sample/Contig/<method...>
    table -- same reshape as Koverage's own coverm_combine rule."""
    input:
        expand(os.path.join(OUTDIR, "preprocess", "temp", "{sample}.cov"), sample=SAMPLE_NAMES)
    output:
        os.path.join(OUTDIR, "preprocess", "results", "sample_coverm_coverage.tsv")
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


rule run_combine_cov:
    """Sample\tContig\tCount\tRPKM\tTPM\tMean\tCovered_fraction\tVariance\n"""
    input:
        os.path.join(OUTDIR, "preprocess", "results", "sample_coverm_coverage.tsv")
    output:
        os.path.join(OUTDIR, "preprocess", "coverage.tsv")
    shell:
        # NR>1 skips the header inside awk itself, rather than the old `sed -i '1d'
        # {input}` which mutated the input in place. That made the rule non-
        # idempotent: a retry after partial failure would see an already-header-
        # stripped input, silently drop its first real data row as if it were still
        # the header, and produce a wrong-but-plausible coverage.tsv rather than
        # erroring. This version never touches {input} at all.
        """
        awk -F '\t' 'NR>1 {{ sum[$2] += $6 }} END {{ for (key in sum) print key, sum[key] }}' {input} > {output}
        """
