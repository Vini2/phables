# Checking the quality of resolved genomes

The sequences of the resolved genomic paths can be found in `resolved_paths.fasta`. Each entry in this FASTA file is a resolved genome (not a contig) and can be directly evaluated using a dedicated viral evaluation tool like [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/).

Here is an example command to run CheckV on the resolved genomes.

```bash
checkv end_to_end resolved_paths.fasta checkv_resolved_paths
```