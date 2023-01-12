"""
Use FragGeneScan and HMMER to scan for bacterial single-copy marker genes in unitigs.
User mmseqs2 to scan for PHROGs in unitigs.
"""

rule scan_smg:
    input:
        {EDGES_FILE}
    params:
        hmm = os.path.join(DBPATH, "marker.hmm")
    output:
        os.path.join(OUTDIR, "edges.fasta.hmmout")
    log:
        frag_out=os.path.join(LOGSDIR, "smg_scan_frag_out.log"),
        frag_err=os.path.join(LOGSDIR, "smg_scan_frag_err.log"),
        hmm_out=os.path.join(LOGSDIR, "smg_scan_hmm_out.log"),
        hmm_err=os.path.join(LOGSDIR, "smg_scan_hmm_err.log")
    conda: 
        "../envs/smg.yaml"
    shell:
        """
            run_FragGeneScan.pl -genome={input} -out={input}.frag -complete=0 -train=complete -thread={threads} 1>{log.frag_out} 2>{log.frag_err}
            hmmsearch --domtblout {output} --cut_tc --cpu {threads} {params.hmm} {input}.frag.faa 1>{log.hmm_out} 2> {log.hmm_err}
        """


rule scan_phrogs:
    input:
        {EDGES_FILE}
    params:
        db = os.path.join(DBPATH, "phrogs_mmseqs_db", "phrogs_profile_db")
    output:
        os.path.join(PHROGS_PATH, "phrogs_annotations.tsv")
    log:
        os.path.join(LOGSDIR, "phrogs_scan.log")
    conda: 
        "../envs/mmseqs.yaml"
    shell:
        """
            mmseqs createdb {input} {PHROGS_PATH}target_seq > {log}
            mmseqs search {PHROGS_PATH}target_seq {params.db} {PHROGS_PATH}results_mmseqs {PHROGS_PATH}tmp/ -s 7 > {log}
            mmseqs createtsv {PHROGS_PATH}target_seq {params.db} {PHROGS_PATH}results_mmseqs {output} --full-header > {log}
        """