# Pre-processing

Phables requires paired-end sequencing reads from viral metagenomic samples to be assembled and preprocessed. The following steps explain the preprocessing steps required to be carried out beforehand.


## Step 1: Assemble the samples

Phables requires the assembly graph file in **Graphical Fragment Assembly (GFA)** format. You can use any assembler that produces the assembly graph in GFA format to assemble your samples OR you can convert a FASTG file to GFA format.

If you have multiple samples, you can do a cross assembly or pool together reads and do a co-assembly.


## Step 2: Run `phables preprocess`

Phables has a `preprocess` subcommand which you can run to process your data. You can see the following command-line options of Phables using `phables preprocess --help`.

```bash
Usage: phables preprocess [OPTIONS] [SNAKE_ARGS]...

  Preprocess data

Options:
  --input PATH                  Path to assembly graph file in .GFA format
                                [required]
  --reads PATH                  Path to directory containing paired-end reads
                                [required]
  --output PATH                 Output directory  [default: phables.out]
  --configfile TEXT             Custom config file [default:
                                (outputDir)/config.yaml]
  --threads INTEGER             Number of threads to use  [default: 1]
  --use-conda / --no-use-conda  Use conda for Snakemake rules  [default: use-
                                conda]
  --conda-prefix PATH           Custom conda env directory
  --snake-default TEXT          Customise Snakemake runtime args  [default:
                                --rerun-incomplete, --printshellcmds,
                                --nolock, --show-failed-logs]
  -h, --help                    Show this message and exit.
```

Assuming your assembly graph file is `assembly_graph.gfa` and reads folder as `fastq`, you can run the `preprocess` subcommand as follows.

```bash
# Preprocess data using 8 threads (default is 1 thread)
phables preprocess --input assembly_graph.gfa --reads fastq --threads 8
```

Note that you should provide the path to the GFA file to the `--input` parameter and the folder containing your sequencing reads to the `--reads` parameter. 

Before running Phables, please make sure that the names of the read files are in the following format. Assuming that your paired-end sequencing reads are in the folder `fastq`, please make sure that the reads are in the format `{sampleName}_R1{fileExtension}`. `fileExtension` can be `.fq`, `.fastq`, `.fq.gz` or `.fastq.gz`.  For example,

```
sample1_R1.fastq.gz
sample1_R2.fastq.gz
sample2_R1.fastq.gz
sample2_R2.fastq.gz
...
```

The output of Phables is set by default to `phables.out`. You can change this using the `--output` argument.

```bash
# Preprocess data using 8 threads (default is 1 thread)
phables preprocess --input assembly_graph.gfa --reads fastq --output my_output_folder --threads 8
```

## Output of preprocessing

The following preprocessing steps will be carried out and your `phables.out` (by default, or `my_output_folder`) output folder will contain the corresponding files and folders.

* Obtain unitig sequences from assembly graph - `edges.fasta`
* Map reads to unitig sequences and get BAM files - `bam_files`
* Run CoverM to get coverage of unitig sequences - `coverage.tsv`
* Scan unitig sequences for single-copy marker genes - `edges.fasta.hmmout`
* Scan unitig sequences for Prokaryotic Virus Remote Homologous Groups (PHROGs) - `phrogs/phrogs_annotations.tsv`

Now we are ready to run Phables.