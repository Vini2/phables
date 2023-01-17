# Phables Usage

## Phables run options

You can see the following command-line options of Phables using `phables run --help`.

```
Usage: phables run [OPTIONS] [SNAKE_ARGS]...

  Run Phables

Options:
  --input PATH                  Path to Hecatomb output  [required]
  --minlength INTEGER           minimum length of circular unitigs to consider
  --mincov INTEGER              minimum coverage of paths to output
  --compcount INTEGER           maximum unitig count to consider a component
  --maxpaths INTEGER            maximum number of paths to resolve for a
                                component
  --mgfrac FLOAT                length threshold to consider single copy
                                marker genes
  --alignscore FLOAT            minimum alignment score for phrog annotations
  --seqidentity FLOAT           minimum sequence identity for phrog
                                annotations
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

  
  DOWNLOAD AND SETUP DATABASES
  phables install
  
  
  PREPROCESSING DATA
  phables preprocess <options> 
  
  
  RUN PHABLES
  phables run <options> 
  For more information on Phables please visit:
  https://phables.readthedocs.io/
  
  
  CLUSTER EXECUTION:
  phables run ... --profile [profile]
  For information on Snakemake profiles see:
  https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
  
  RUN EXAMPLES:
  Required:           phables run --input [folder]
  Specify threads:    phables run ... --threads [threads]
  Disable conda:      phables run ... --no-use-conda 
  Change defaults:    phables run ... --snake-default="-k --nolock"
  Add Snakemake args: phables run ... --dry-run --keep-going --touch
  Specify targets:    phables run ... all print_targets
  Available targets:
      all             Run everything (default)
      print_targets   List available targets
```

## Example usage

```bash
phables run --input hecatomb.out
```

Note that you should provide the path to the Hecatomb output folder `hecatomb.out` for `--input`.

The output of all the Phables subcommands is set by default to `phables.out`. If you changed the `--output` path in the [preprocessing steps](https://phables.readthedocs.io/en/latest/preprocess/) to `my_output_folder` for example, make sure to change the `--output` parameter for `phables run` as follows.

```bash
phables run --input hecatomb.out  --output my_output_folder
```

## Output

The output of Phables will contain the following main files and folders.

* `resolved_paths.fasta` containing the resolved genomes
* `resolved_phages` folder containing the resolved genomes in individual FASTA files
* `resolved_genome_info.txt` containing the path name, coverage, length, GC content and unitig order of the resolved genomes
* `resolved_edges.fasta` containing the unitigs that make up the resolved genomes
* `resolved_component_info.txt` containing the details of the phage bubbles resolved

## Reporting Issues

Phables is still under testing. If you want to test (or break) Phables give it a try and report any issues and suggestions under [Phables Issues](https://github.com/Vini2/phables/issues).