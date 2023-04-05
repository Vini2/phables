"""
Declare your targets here!
A separate file is ideal if you have lots of target files to create, or need some python logic to determine
the targets to declare. This example shows targets that are dependent on the input file type.
"""

allTargets = []

allTargets.append(os.path.join(TESTDIR, "resolved_paths.fasta"))
allTargets.append(os.path.join(TESTDIR, "resolved_genome_info.txt"))
allTargets.append(os.path.join(TESTDIR, "resolved_edges.fasta"))
allTargets.append(os.path.join(TESTDIR, "resolved_component_info.txt"))
allTargets.append(os.path.join(TESTDIR, "phage_like_edges.fasta"))