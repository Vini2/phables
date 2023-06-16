"""
Use Minimap2 to map the reads of all the samples to unitigs in the edges.fasta 
file from step 2 and Samtools to index the BAM files.
"""

# rule map_reads_to_unitigs:
#     input:
#         edges = EDGES_FILE,
#         r1 = os.path.join(READ_DIR, PATTERN_R1),
#         r2 = os.path.join(READ_DIR, PATTERN_R2)
#     output:
#         bam = os.path.join(BAM_PATH, "{sample}.bam"),
#         index = os.path.join(BAM_PATH, "{sample}.bam.bai")
#     threads:
#         config["resources"]["jobCPU"]
#     resources:
#         mem_mb = config["resources"]["jobMem"]
#     log:
#         os.path.join(LOGSDIR, "{sample}_mapping.log")
#     conda:
#         os.path.join("..", "envs", "mapping.yaml")
#     shell:
#         """
#             minimap2 -t {threads} -N 5 -ax sr {input.edges} {input.r1} {input.r2} | samtools sort -@ {threads} > {output.bam}
#             samtools index -@ {threads} {output.bam} {output.index}
#         """