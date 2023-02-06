
preprocessTargets = []
phablesTargets = []


"""PREPROCESSING TARGETS"""
EDGES_FILE = os.path.join(OUTDIR, "edges.fasta")
preprocessTargets.append(EDGES_FILE)

BAM_PATH = os.path.join(OUTDIR, 'bam_files/')
preprocessTargets.append(expand(os.path.join(BAM_PATH, "{sample}.bam"), sample=SAMPLES))
preprocessTargets.append(expand(os.path.join(BAM_PATH, "{sample}.bam.bai"), sample=SAMPLES))

COVERM_PATH = os.path.join(OUTDIR, 'coverage_rpkm/')
preprocessTargets.append(expand(os.path.join(COVERM_PATH, "{sample}_rpkm.tsv"), sample=SAMPLES))
preprocessTargets.append(os.path.join(OUTDIR, "coverage.tsv"))
preprocessTargets.append(os.path.join(OUTDIR, "edges.fasta.hmmout"))

PHROGS_PATH = os.path.join(OUTDIR, 'phrogs/')
preprocessTargets.append(os.path.join(PHROGS_PATH, "phrogs_annotations.tsv"))


"""PHABLES TARGETS"""
phablesTargets.append(os.path.join(OUTDIR, "resolved_paths.fasta"))
phablesTargets.append(os.path.join(OUTDIR, "resolved_genome_info.txt"))
phablesTargets.append(os.path.join(OUTDIR, "resolved_edges.fasta"))
phablesTargets.append(os.path.join(OUTDIR, "resolved_component_info.txt"))
