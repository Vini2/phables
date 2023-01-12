# Pre-processing

Phables requires the sequencing reads from viral metagenomic samples to be run through [Hecatomb](https://hecatomb.readthedocs.io/en/latest/) and be preprocessed. The following steps explain the preprocessing steps required to be carried out before.


## Step 1: Run read samples through Hecatomb

Phables requires the assembly output from [Hecatomb](https://hecatomb.readthedocs.io/en/latest/). First, run your sequencing read samples through Hecatomb using the following command.

```bash
hecatomb run --reads <reads_folder>
```

You can find the assembly output in `hecatomb.out/processing/assembly/CONTIG_DICTIONARY/FLYE/`. Phables will be using the `assembly_graph.gfa` and `assembly_info.txt` files.

## Step 2: Run `phables preprocess`

Phables has a `preprocess` subcommand which you can run to process your data as follows.

```bash
phables preprocess --input <path_to_hecatomb.out> --reads <reads_folder>
```

Note that you should provide the path to the Hecatomb output folder `hecatomb.out` for `--input` and the reads folder to `--reads`.

The output of Phables is set by default to `phables.out`. You can change this using the `--output` argument.

The following preprocessing steps will be carried out by the `preprocess` commands.

* Obtain unitig sequences from assembly graph
* Map reads to unitig sequences and get BAM files
* Run CoverM to get coverage of unitig sequences
* Scan unitig sequences for single-copy marker genes and PHROGs

Now we are ready to run Phables.