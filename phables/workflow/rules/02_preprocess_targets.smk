"""
Declare your targets here!
A separate file is ideal if you have lots of target files to create, or need some python logic to determine
the targets to declare. This example shows targets that are dependent on the input file type.
"""

preprocess_files = []

preprocess_files.append(os.path.join(OUTDIR, "edges.fasta"))