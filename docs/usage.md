# Phables Usage

Phables run options can be found using the `phables run -h` command.

```
Usage: phables run [OPTIONS] [SNAKE_ARGS]...

  Run Phables

Options:
  --input PATH                  Path to assembly graph file in .GFA format
                                [required]
  --reads PATH                  Path to directory containing paired-end reads
                                [required]
  --minlength INTEGER           minimum length of circular unitigs to consider
                                [default: 2000]
  --mincov INTEGER              minimum coverage of paths to output  [default:
                                10]
  --compcount INTEGER           maximum unitig count to consider a component
                                [default: 200]
  --maxpaths INTEGER            maximum number of paths to resolve for a
                                component  [default: 10]
  --mgfrac FLOAT                length threshold to consider single copy
                                marker genes  [default: 0.2]
  --evalue FLOAT                maximum e-value for phrog annotations
                                [default: 1e-10]
  --seqidentity FLOAT           minimum sequence identity for phrog
                                annotations  [default: 0.3]
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

## Run options explained

* `--input` - assembly graph file in .GFA format
* `--reads` - folder containing paired-end read files
* `--minlength` - minimum length of circular unitigs to consider [default: 2000]
* `--mincov` - minimum coverage of paths to output [default: 10]
* `--compcount` - maximum unitig count to consider a component [default: 200]
* `--maxpaths` - maximum number of paths to resolve for a component [default: 10]
* `--mgfrac` - length threshold to consider single copy marker genes [default: 0.2]
* `--evalue` - maximum e-value for phrog annotations [default: 1e-10]
* `--seqidentity` - minimum sequence identity for phrog annotations [default: 0.3]
* `--output` - path to the output directory [default: `phables.out`]
* `--configfile` - custom config file [default: `(outputDir)/config.yaml`]
* `--threads` - number of threads to use  [default: 1]
* `--use-conda` / `--no-use-conda` - use conda for Snakemake rules  [default: `use-conda`]
* `--conda-prefix` - custom conda env directory
* `--snake-default` - customise Snakemake runtime args  [default: `--rerun-incomplete, --printshellcmds, --nolock, --show-failed-logs`]


## Example usage

Assuming your assembly graph file is `assembly_graph.gfa` and reads folder as `fastq`, you can run `phables` as follows.

```bash
# Preprocess data using 8 threads (default is 1 thread)
phables run --input assembly_graph.gfa --reads fastq --threads 8
```

Note that you should provide the path to the GFA file to the `--input` parameter and the folder containing your sequencing reads to the `--reads` parameter. 

The output of Phables is set by default to `phables.out`. You can update the output path using the `--output` parameter for `phables run` as follows.

```bash
# Preprocess data using 8 threads (default is 1 thread)
phables run --input assembly_graph.gfa --reads fastq  --output my_output_folder --threads 8
```

The `phables run` command will run the following preprocessing steps with the corresponding files and folders, prior to genome resolution.

* Obtain unitig sequences from assembly graph - `edges.fasta`
* Map reads to unitig sequences and get BAM files - `bam_files`
* Run CoverM to get coverage of unitig sequences - `coverage.tsv`
* Scan unitig sequences for single-copy marker genes - `edges.fasta.hmmout`
* Scan unitig sequences for Prokaryotic Virus Remote Homologous Groups (PHROGs) - `phrogs/phrogs_annotations.tsv`

You can use the following command to run only the above preprocessing steps.

```
phables run preprocessing --input assembly_graph.gfa --reads fastq --threads 8
```

## Output

The output of Phables will contain the following main files and folders.

* `resolved_paths.fasta` containing the resolved genomes
* `resolved_phages` folder containing the resolved genomes in individual FASTA files
* `resolved_genome_info.txt` containing the path name, coverage, length, GC content and unitig order of the resolved genomes
* `resolved_edges.fasta` containing the unitigs that make up the resolved genomes
* `resolved_component_info.txt` containing the details of the phage bubbles resolved