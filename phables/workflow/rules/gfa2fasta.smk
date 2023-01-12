"""
Run gfa2fasta to obtain the sequences corresponding to edges in the Flye and Miniasm assembly graphs in FASTA format.
The assembly graph file of Flye (assembly_graph.gfa) should be provided as inputs.
"""

rule run_gfa2fasta:
    input:
        GRAPH_FILE
    output:
        EDGES_FILE
    params:
        assembler = "flye",
        graph = GRAPH_FILE,
        output = OUTDIR
    log:
        os.path.join(LOGSDIR, "gfa2fasta.log")
    conda: 
        "../envs/phables.yaml"
    script:
        os.path.join('..', 'scripts', 'gfa2fasta.py')