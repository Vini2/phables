# Phables Usage

## Phables options

You can see the following command-line options of Phables using `python phables --help`.

```
Usage: phables [OPTIONS]

  Phables: Phage bubbles resolve bacteriophage genomes in viral metagenomic
  samples. Please refer the full documentation available on Read the Docs at
  https://phables.readthedocs.io/

Options:
  -g, --graph PATH          path to the assembly graph file  [required]
  -p, --paths PATH          path to the assembly path info file  [required]
  -c, --coverage PATH       path to the coverage file  [required]
  -b, --bampath PATH        path to the bam files  [required]
  -hm, --hmmout PATH        path to the .hmmout file  [required]
  -ph, --phrogs PATH        path to the phrog annotations file  [required]
  -ml, --minlength INTEGER  minimum length of circular unitigs to consider
  -mcov, --mincov INTEGER   minimum coverage of paths to output
  -cc, --compcount INTEGER  maximum unitig count to consider a component
  -mp, --maxpaths INTEGER   maximum number of paths to resolve for a component
  -mgf, --mgfrac FLOAT      length threshold to consider single copy marker
                            genes
  -as, --alignscore FLOAT   minimum alignment score for phrog annotations
  -si, --seqidentity FLOAT  minimum sequence identity for phrog annotations
  -o, --output PATH         path to the output folder  [required]
  --version                 Show the version and exit.
  --help                    Show this message and exit.
```

## Input

Phables takes the following files/paths as input.

* Assembly graph file (`assembly_graph.gfa`)
* Path information of contigs (`assembly_info.txt`)
* Path to BAM files (`bam_files/`) from preprocessing step 3
* Coverage of unitig sequences in each sample in `.tsv` format (`coverage.tsv`) from preprocessing step 4
* HMMER result for scan of single-copy marker genes (`edges.fasta.hmmout`) from preprocessing step 5
* MMseqs result for scan of PHROGs (`phrog_annot.tsv`) from preprocessing step 5

## Example usage

```bash
phables -g assembly_graph.gfa -p assembly_info.txt -hm edges.fasta.hmmout -ph phrog_annot.tsv -c coverage.tsv -b bam_files/ -o /output/path/
```

## Output

The output from Phables will be as follows.

* `resolved_paths.fasta` containing the resolved genomes
* `resolved_phages` folder containing the resolved genomes in individual FASTA files
* `resolved_genome_info.txt` containing the path name, coverage, length, GC content and unitig order of the resolved genomes
* `resolved_edges.fasta` containing the unitigs that make up the resolved genomes
* `resolved_component_info.txt` containing the details of the phage bubbles resolved

## Snakemake pipeline
To Do

## Reporting Issues

Phables is still under testing. If you want to test (or break) Phables give it a try and report any issues and suggestions under [Phables Issues](https://github.com/Vini2/phables/issues).