"""
Run gfa2fasta to obtain the sequences corresponding to unitigs in the assembly graphs in FASTA format.
The assembly graph file with .GFA extension should be provided as inputs.
"""

rule run_gfa2fasta:
    input:
        GRAPH_FILE
    output:
        EDGES_FILE
    params:
        graph = GRAPH_FILE,
        output = os.path.join(OUTDIR, "preprocess"),
        log = os.path.join(LOGSDIR, "gfa2fasta.log")
    log:
        os.path.join(LOGSDIR, "gfa2fasta.log")
    conda: 
        os.path.join("..", "envs", "phables.yaml")
    script:
        os.path.join('..', 'scripts', 'gfa2fasta.py')