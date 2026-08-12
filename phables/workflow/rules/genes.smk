"""
Call genes on unitigs, then use HMMER to scan for bacterial single-copy marker genes.
Use mmseqs2 to scan for PHROGs in unitigs.

Gene calling is a separate rule from the marker-gene search so the caller can be
swapped via --genecaller without touching anything downstream. Both callers emit
the same `{seqid}_{start}_{end}_{strand}` protein id convention, which
gene_utils.get_smg_unitigs relies on to recover the unitig name.
"""

PROTEINS_FILE = EDGES_FILE + ".frag.faa"


if GC == "pyrodigal-gv":

    rule call_genes:
        input:
            genome = EDGES_FILE,
        threads:
            config["resources"]["jobCPU"]
        resources:
            mem_mb = config["resources"]["jobMem"]
        output:
            faa = PROTEINS_FILE
        log:
            os.path.join(LOGSDIR, "gene_call_pyrodigal_gv.log")
        conda:
            None if CONTAINER_IMAGE else os.path.join("..", "envs", "genecall.yaml")
        container:
            CONTAINER_IMAGE
        script:
            os.path.join("..", "scripts", "gene_caller.py")

else:

    rule call_genes:
        input:
            genome = EDGES_FILE,
        threads:
            config["resources"]["jobCPU"]
        resources:
            mem_mb = config["resources"]["jobMem"]
        output:
            faa = PROTEINS_FILE
        params:
            frag = EDGES_FILE + ".frag",
        log:
            out = os.path.join(LOGSDIR, "gene_call_fraggenescan_out.log"),
            err = os.path.join(LOGSDIR, "gene_call_fraggenescan_err.log"),
        conda:
            None if CONTAINER_IMAGE else os.path.join("..", "envs", "smg.yaml")
        container:
            CONTAINER_IMAGE
        shell:
            """
                run_FragGeneScan.pl -genome={input.genome} -out={params.frag} -complete=0 -train=complete -thread={threads} 1>{log.out} 2>{log.err}
            """


rule scan_smg:
    input:
        faa = PROTEINS_FILE,
        hmm = os.path.join(DBPATH, "marker.hmm"),
    threads:
        config["resources"]["jobCPU"]
    resources:
        mem_mb = config["resources"]["jobMem"]
    output:
        hmmout = os.path.join(OUTDIR, "preprocess", "edges.fasta.hmmout")
    log:
        hmm_out=os.path.join(LOGSDIR, "smg_scan_hmm_out.log"),
        hmm_err=os.path.join(LOGSDIR, "smg_scan_hmm_err.log")
    conda:
        None if CONTAINER_IMAGE else os.path.join("..", "envs", "smg.yaml")
    container:
        CONTAINER_IMAGE
    shell:
        """
            hmmsearch --domtblout {output.hmmout} --cut_tc --cpu {threads} {input.hmm} {input.faa} 1>{log.hmm_out} 2> {log.hmm_err}
        """


rule scan_phrogs:
    input:
        genome = EDGES_FILE,
        db = os.path.join(DBPATH,"phrogs_mmseqs_db","phrogs_profile_db")
    threads:
        config["resources"]["jobCPU"]
    resources:
        mem_mb = config["resources"]["jobMem"]
    output:
        os.path.join(OUTDIR, "preprocess", "phrogs_annotations.tsv")
    params:
        out_path = os.path.join(OUTDIR, "preprocess", "phrogs"),
        target_seq = os.path.join(OUTDIR, "preprocess", "phrogs", "target_seq"),
        results_mmseqs = os.path.join(OUTDIR, "preprocess", "phrogs", "results_mmseqs"),
        tmp = os.path.join(OUTDIR, "preprocess", "phrogs", "tmp"),
    log:
        os.path.join(LOGSDIR, "phrogs_scan.log")
    conda:
        None if CONTAINER_IMAGE else os.path.join("..", "envs", "mmseqs.yaml")
    container:
        CONTAINER_IMAGE
    shell:
        """
        mkdir -p {params.out_path}
        mmseqs createdb {input} {params.target_seq} > {log}
        mmseqs search {params.target_seq} {input.db} {params.results_mmseqs} {params.tmp} --threads {threads} -s 7 > {log}
        mmseqs createtsv {params.target_seq} {input.db} {params.results_mmseqs} {output} --threads {threads} --full-header > {log}
        rm -rf {params.out_path}
        """


"""
Structural phage-gene detection: ProstT5 -> 3Di -> foldseek search against a
hallmark PHROG structure subDB, selected via --phage-detection prostt5-foldseek.

Gated behind `if PD == "prostt5-foldseek"`, unlike scan_phrogs above -- not just an
optimisation. hallmark_db/prostt5_checkpoint/etc default to empty in config.yaml
when unused, and Snakemake parses every rule's input/output at load time regardless
of whether it ends up in the target list, so an un-gated `hallmark_db = config["hallmark_db"]`
input fails with `SyntaxError: Input and output files have to be specified as
strings or lists of strings` on a plain `phables run` even when this path is never
selected. Only relying on the target list (as scan_phrogs does) isn't enough here
because these params, not just outputs, can be structurally invalid at parse time.

These consume PROTEINS_FILE, the same protein set the marker-gene search uses
(genes.smk top), unifying what used to be two independent detection paths (see
the fork audit, section 2.1: scan_phrogs previously ran mmseqs directly on the
nucleotide edges FASTA and never shared gene calls with the marker-gene path).
"""

if PD == "prostt5-foldseek":

    QUERY_3DI = os.path.join(OUTDIR, "preprocess", "hallmark", "proteins_3di.fasta")


    if PROSTT5_CONTAINER:

        # container:-only, deliberately with no conda: alongside it. Snakemake's
        # documented way to combine the two ("Ad-hoc combination of Conda package
        # management with containers") builds a fresh, isolated conda env *inside*
        # the container rather than exposing what the image already has installed
        # -- which would just reinstall a second, redundant torch and ignore the
        # image's own verified-working one (e.g. phold's own container, which
        # already bundles pholdlib + a working torch/ROCm stack for Setonix, per
        # the same logic phold's own Snakemake rules use). Needs
        # `--use-singularity` passed through phables' snake_args passthrough --
        # --use-conda alone won't honour this directive.
        rule predict_3di:
            input:
                faa = PROTEINS_FILE,
            threads:
                config["resources"]["jobCPU"]
            resources:
                mem_mb = config["resources"]["jobMem"]
            output:
                threedi = QUERY_3DI
            params:
                checkpoint = config["prostt5_checkpoint"],
                model_name = config["prostt5_model"],
                model_dir = config["prostt5_model_dir"],
                half_precision = config["prostt5_half_precision"],
                cpu = config["prostt5_cpu"],
                max_residues = config["prostt5_max_residues"],
                max_seq_len = config["prostt5_max_seq_len"],
                max_batch = config["prostt5_max_batch"],
            log:
                os.path.join(LOGSDIR, "predict_3di.log")
            container:
                PROSTT5_CONTAINER
            script:
                os.path.join("..", "scripts", "predict_3di.py")

    else:

        rule predict_3di:
            input:
                faa = PROTEINS_FILE,
            threads:
                config["resources"]["jobCPU"]
            resources:
                mem_mb = config["resources"]["jobMem"]
            output:
                threedi = QUERY_3DI
            params:
                checkpoint = config["prostt5_checkpoint"],
                model_name = config["prostt5_model"],
                model_dir = config["prostt5_model_dir"],
                half_precision = config["prostt5_half_precision"],
                cpu = config["prostt5_cpu"],
                max_residues = config["prostt5_max_residues"],
                max_seq_len = config["prostt5_max_seq_len"],
                max_batch = config["prostt5_max_batch"],
            log:
                os.path.join(LOGSDIR, "predict_3di.log")
            conda:
                # gpu_backend selects which torch build this env solves against --
                # cpu/cuda/rocm need different PyTorch wheels (conda envs are
                # solved once from a static file, so this has to be three files,
                # not one file with a runtime switch). See the individual env
                # files for what each backend actually needs and why.
                os.path.join("..", "envs", f"prostt5-{GPU_BACKEND}.yaml")
            script:
                os.path.join("..", "scripts", "predict_3di.py")


    rule build_hallmark_query_db:
        input:
            faa = PROTEINS_FILE,
            threedi = QUERY_3DI,
        output:
            db = os.path.join(OUTDIR, "preprocess", "hallmark", "query_db"),
        params:
            out_prefix = os.path.join(OUTDIR, "preprocess", "hallmark", "query_db"),
        log:
            os.path.join(LOGSDIR, "build_hallmark_query_db.log")
        conda:
            None if CONTAINER_IMAGE else os.path.join("..", "envs", "foldseek.yaml")
        container:
            CONTAINER_IMAGE
        script:
            os.path.join("..", "scripts", "build_foldseek_query_db.py")


    if FOLDSEEK_GPU and GPU_BACKEND != "cuda":
        raise ValueError(
            "foldseek_gpu is True but gpu_backend is "
            f"'{GPU_BACKEND}', not 'cuda'. foldseek's --gpu mode is CUDA-only "
            "(confirmed against foldseek's own README/source -- there is no ROCm "
            "or Metal build); it cannot accelerate on an AMD (ROCm, e.g. Setonix's "
            "MI250X) or CPU backend. Set gpu_backend: cuda, or turn foldseek_gpu off."
        )

    # Bioconda's foldseek is CPU-only. GPU search needs foldseek's separate
    # foldseek-linux-gpu.tar.gz build (not a conda package at all -- see
    # foldseek's README "Installation" section) and a target DB reformatted
    # with `foldseek makepaddedseqdb`, conventionally named with a `_gpu`
    # suffix -- both details taken directly from phold's own GPU search code
    # (phold/features/run_foldseek.py), which this mirrors rather than
    # reimplementing from the foldseek README alone (the README suggests one
    # padded DB works for both CPU and GPU; phold's actual working code
    # maintains a separate _gpu-suffixed DB, which is the safer bet to follow
    # here since it's the reference implementation in this ecosystem).
    #
    # NOT independently verified against real CUDA hardware -- this machine
    # has neither an NVIDIA GPU nor Setonix access. The plumbing (flags, DB
    # suffix, env selection) is wired correctly per phold's own precedent;
    # whether it actually accelerates anything is unconfirmed.
    HALLMARK_TARGET_DB = (
        f"{config['hallmark_db']}_gpu" if FOLDSEEK_GPU else config["hallmark_db"]
    )
    _foldseek_gpu_flags = "--gpu 1 --prefilter-mode 1" if FOLDSEEK_GPU else ""


    rule scan_hallmark:
        input:
            db = os.path.join(OUTDIR, "preprocess", "hallmark", "query_db"),
            hallmark_db = HALLMARK_TARGET_DB,
        threads:
            config["resources"]["jobCPU"]
        resources:
            mem_mb = config["resources"]["jobMem"]
        output:
            os.path.join(OUTDIR, "preprocess", "hallmark_hits.tsv")
        params:
            query_prefix = os.path.join(OUTDIR, "preprocess", "hallmark", "query_db"),
            result = os.path.join(OUTDIR, "preprocess", "hallmark", "result"),
            tmp = os.path.join(OUTDIR, "preprocess", "hallmark", "tmp"),
            gpu_flags = _foldseek_gpu_flags,
        log:
            os.path.join(LOGSDIR, "scan_hallmark.log")
        conda:
            None if CONTAINER_IMAGE else os.path.join("..", "envs", "foldseek.yaml")
        container:
            CONTAINER_IMAGE
        shell:
            """
            foldseek search {params.query_prefix} {input.hallmark_db} {params.result} {params.tmp} \
                --threads {threads} -s 7 -e {config[hallmark_evalue]} {params.gpu_flags} > {log}
            foldseek convertalis {params.query_prefix} {input.hallmark_db} {params.result} {output} \
                --format-output query,target,evalue,bits,fident --threads {threads} >> {log}
            rm -rf {params.tmp}
            """