# Pre-processing

Phables requires paired-end sequencing reads from viral metagenomic samples to be run through [Hecatomb](https://hecatomb.readthedocs.io/en/latest/) and be preprocessed. The following steps explain the preprocessing steps required to be carried out beforehand.


## Step 1: Run read samples through Hecatomb

Phables requires the assembly output from [Hecatomb](https://hecatomb.readthedocs.io/en/latest/). First, you have to run your sequencing read samples through Hecatomb (version 1.1.0 or higher). 

Assume that your paired-end sequencing reads are in the folder `fastq`. Make sure that the reads are in the format `{sampleName}_R1{fileExtension}`. `fileExtension` can be `.fq`, `.fastq`, `.fq.gz` or `.fastq.gz`.  For example,

```
sample1_R1.fastq.gz
sample1_R2.fastq.gz
sample2_R1.fastq.gz
sample2_R2.fastq.gz
...
```

Now we can run Hecatomb using the following command.

```bash
hecatomb run --reads fastq/
```

The hecatomb output will be saved in `hecatomb.out`.

You can find the assembly output in `hecatomb.out/processing/assembly/CONTIG_DICTIONARY/FLYE/` and Phables will be using the `assembly_graph.gfa` and `assembly_info.txt` files found in this path.

## Step 2: Run `phables preprocess`

Phables has a `preprocess` subcommand which you can run to process your data. Assuming your hecatomb output folder as `hecatomb.out` and reads folder as `fastq`, you can run the `preprocess` subcommand as follows.

```bash
# Preprocess data using 8 threads (default is 1 thread)
phables preprocess --input hecatomb.out --reads fastq --threads 8
```

Note that you should provide the path to the Hecatomb output folder `hecatomb.out` to the `--input` parameter and the folder containing your sequencing reads to the `--reads` parameter.

The output of Phables is set by default to `phables.out`. You can change this using the `--output` argument.

```bash
# Preprocess data using 8 threads (default is 1 thread)
phables preprocess --input hecatomb.out --reads fastq --output my_output_folder --threads 8
```

The following preprocessing steps will be carried out by the `preprocess` commands.

* Obtain unitig sequences from assembly graph
* Map reads to unitig sequences and get BAM files
* Run CoverM to get coverage of unitig sequences
* Scan unitig sequences for single-copy marker genes and PHROGs

## Output

If you ran Phables with the default output path, ypur `phables.out` output folder will have the following files.



Now we are ready to run Phables.