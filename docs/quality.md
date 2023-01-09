# Checking the quality of resolved genomes

The sequences of the resolved genomic paths can be found in `resolved_paths.fasta`. Each entry in this FASTA file is a resolved genome (not a contig) and can be directly evaluated using a dedicated viral evaluation tool like [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/). The following sections will walk you through how to setup and run CheckV.

## Installing CheckV

The recommended way to install CheckV is using [`conda`](https://docs.conda.io/en/latest/).

```bash
# Create a new conda environment and install checkv
conda create -n checkv -c conda-forge -c bioconda checkv

# Activate checkv conda environment
conda activate checkv
```

You can also install using [`pip`](https://pip.pypa.io/en/stable/).

```bash
pip install checkv
```

## Download the CheckV database

```bash
checkv download_database ./
```

Now you need to to specify the `CHECKVDB` location.

```bash
export CHECKVDB=/path/to/checkv-db
```

## Running CheckV

Here is an example command to run CheckV on the resolved genomes.

```bash
checkv end_to_end resolved_paths.fasta checkv_resolved_paths -t 16
```

The `end_to_end` option will run the full pipeline.

You can also run individual commands for each step in the pipeline as follows.

```bash
checkv contamination resolved_paths.fasta checkv_resolved_paths -t 16
checkv completeness resolved_paths.fasta checkv_resolved_paths -t 16
checkv complete_genomes resolved_paths.fasta checkv_resolved_paths
checkv quality_summary resolved_paths.fasta checkv_resolved_paths
```

## CheckV outputs

CheckV will produce the following `.tsv` files.

* `complete_genomes.tsv` - overview of putative complete genomes identified
* `completeness.tsv` - overview of how completeness was estimated
* `contamination.tsv` - overview of how contamination was estimated
* `quality_summary.tsv` - integrated quality results