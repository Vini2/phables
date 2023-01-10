"""
gfa2fasta.py: Obtain the sequences corresponding to edges in the Flye and Miniasm assembly graphs in FASTA format.
The assembly graph file of Flye (assembly_graph.gfa) should be provided as inputs.
"""

rule gfa2fasta:
    input:
        GRAPH_FILE
    output:
        os.path.join(OUTDIR, "edges.fasta"),
    params:
        assembler = "flye"
    log:
        os.path.join(OUTDIR, "gfa2fasta.log")
    conda: 
        "../envs/phables.yaml"
    threads: 1
    shell:
        """
            gfa2fasta --graph {input} --assembler {params.assembler} --output {OUTDIR}
        """