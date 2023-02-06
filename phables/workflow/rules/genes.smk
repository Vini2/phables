"""
Use FragGeneScan and HMMER to scan for bacterial single-copy marker genes in unitigs.
User mmseqs2 to scan for PHROGs in unitigs.
"""

rule scan_smg:
    input:
        genome = EDGES_FILE,
        hmm = os.path.join(DBPATH, "marker.hmm"),
    threads:
        THREADS
    output:
        frag = EDGES_FILE + ".frag",
        frag_faa = EDGES_FILE + ".frag.faa",
        hmmout = os.path.join(OUTDIR, "edges.fasta.hmmout")
    log:
        frag_out=os.path.join(LOGSDIR, "smg_scan_frag_out.log"),
        frag_err=os.path.join(LOGSDIR, "smg_scan_frag_err.log"),
        hmm_out=os.path.join(LOGSDIR, "smg_scan_hmm_out.log"),
        hmm_err=os.path.join(LOGSDIR, "smg_scan_hmm_err.log")
    conda: 
        os.path.join("..", "envs", "smg.yaml")
    shell:
        """
            run_FragGeneScan.pl -genome={input.genome} -out={output.frag} -complete=0 -train=complete -thread={threads} 1>{log.frag_out} 2>{log.frag_err}
            hmmsearch --domtblout {output.hmmout} --cut_tc --cpu {threads} {input.hmm} {output.frag_faa} 1>{log.hmm_out} 2> {log.hmm_err}
        """


rule scan_phrogs:
    input:
        genome = EDGES_FILE,
        db = os.path.join(DBPATH,"phrogs_mmseqs_db","phrogs_profile_db")
    threads:
        THREADS
    output:
        os.path.join(PHROGS_PATH, "phrogs_annotations.tsv")
    params:
        target_seq = os.path.join(PHROGS_PATH, "target_seq"),
        results_mmseqs = os.path.join(PHROGS_PATH, "results_mmseqs"),
        tmp = os.path.join(PHROGS_PATH, "tmp"),
    log:
        os.path.join(LOGSDIR, "phrogs_scan.log")
    conda: 
        os.path.join("..", "envs", "mmseqs.yaml")
    shell:
        """
            mmseqs createdb {input} {params.target_seq} > {log}
            mmseqs search {params.target_seq} {input.db} {params.results_mmseqs} {params.tmp} --threads {threads} -s 7 > {log}
            mmseqs createtsv {params.target_seq} {input.db} {params.results_mmseqs} {output} --threads {threads} --full-header > {log}
        """