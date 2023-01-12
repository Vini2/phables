"""
Declare your targets here!
A separate file is ideal if you have lots of target files to create, or need some python logic to determine
the targets to declare. This example shows targets that are dependent on the input file type.
"""

allTargets = []

EDGES_FILE = os.path.join(OUTDIR, "edges.fasta")
allTargets.append(EDGES_FILE)

BAM_PATH = os.path.join(OUTDIR, 'bam_files/')
allTargets.append(expand(os.path.join(BAM_PATH, "{sample}.bam"), sample=SAMPLES))
allTargets.append(expand(os.path.join(BAM_PATH, "{sample}.bam.bai"), sample=SAMPLES))

COVERM_PATH = os.path.join(OUTDIR, 'coverage_rpkm/')
allTargets.append(expand(os.path.join(COVERM_PATH, "{sample}_rpkm.tsv"), sample=SAMPLES))

allTargets.append(os.path.join(OUTDIR, "coverage.tsv"))

allTargets.append(os.path.join(OUTDIR, "edges.fasta.hmmout"))

PHROGS_PATH = os.path.join(OUTDIR, 'phrogs/')
allTargets.append(os.path.join(PHROGS_PATH, "phrogs_annotations.tsv"))