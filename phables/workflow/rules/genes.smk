"""
Use FragGeneScan and HMMER to scan for bacterial single-copy marker genes in unitigs.
User mmseqs2 to scan for PHROGs in unitigs.
"""

rule scan_smg:
    input:
        genome = EDGES_FILE,
        hmm = os.path.join(DBPATH, "marker.hmm"),
    threads:
        config["resources"]["jobCPU"]
    resources:
        mem_mb = config["resources"]["jobMem"],
        mem = str(config["resources"]["jobMem"]) + "MB"
    output:
        hmmout = os.path.join(OUTDIR, "preprocess", "edges.fasta.hmmout")
    params:
        frag = EDGES_FILE + ".frag",
        frag_faa = EDGES_FILE + ".frag.faa",
    log:
        frag_out=os.path.join(LOGSDIR, "smg_scan_frag_out.log"),
        frag_err=os.path.join(LOGSDIR, "smg_scan_frag_err.log"),
        hmm_out=os.path.join(LOGSDIR, "smg_scan_hmm_out.log"),
        hmm_err=os.path.join(LOGSDIR, "smg_scan_hmm_err.log")
    conda: 
        os.path.join("..", "envs", "smg.yaml")
    shell:
        """
            run_FragGeneScan.pl -genome={input.genome} -out={params.frag} -complete=0 -train=complete -thread={threads} 1>{log.frag_out} 2>{log.frag_err}
            hmmsearch --domtblout {output.hmmout} --cut_tc --cpu {threads} {input.hmm} {params.frag_faa} 1>{log.hmm_out} 2> {log.hmm_err}
        """


rule scan_phrogs:
    input:
        genome = EDGES_FILE,
        db = os.path.join(DBPATH,"phrogs_mmseqs_db","phrogs_profile_db")
    threads:
        config["resources"]["jobCPU"]
    resources:
        mem_mb = config["resources"]["jobMem"],
        mem = str(config["resources"]["jobMem"]) + "MB"
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
        os.path.join("..", "envs", "mmseqs.yaml")
    shell:
        """
        mkdir -p {params.out_path}
        mmseqs createdb {input} {params.target_seq} > {log}
        mmseqs search {params.target_seq} {input.db} {params.results_mmseqs} {params.tmp} --threads {threads} -s 7 > {log}
        mmseqs createtsv {params.target_seq} {input.db} {params.results_mmseqs} {output} --threads {threads} --full-header > {log}
        rm -rf {params.out_path}
        """