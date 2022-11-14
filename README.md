<p align="center">
  <img src="phables_logo.png" width="600" title="phables logo" alt="phables logo">
</p>

Phables: Phage bubbles resolve bacteriophage genomes in viral metagenomic samples
===============

[![CI](https://github.com/Vini2/phables/actions/workflows/testing.yml/badge.svg)](https://github.com/Vini2/phables/actions/workflows/testing.yml)
![GitHub](https://img.shields.io/github/license/Vini2/phables)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Phables is a tool developed to resolve bacteriophage genomes using phage bubbles in viral metagenomic data. It models cyclic phage-like components in the viral metagenomic assembly as flow networks, models as a minimum flow decomposition problem and resolves genomic paths corresponding to flow paths determined. Phables uses the [Minimum Flow Decomposition via  Integer Linear Programming](https://github.com/algbio/MFD-ILP) implementation to obtain the flow paths.

![](Phables_workflow.png)

Table of Contents
-----------
- [Setting up Phables](#setting-up-phables)
  - [Downloading Phables](#downloading-phables)
  - [Installing Phables](#installing-phables)
  - [Setting up Gurobi](#setting-up-gurobi)
  - [Test the setup](#test-the-setup)
- [Pre-processing](#pre-processing)
  - [Step 1: Run read samples through Hecatomb](#step-1-run-read-samples-through-hecatomb)
  - [Step 2: Obtain unitig sequences from assembly graph](#step-2-obtain-unitig-sequences-from-assembly-graph)
  - [Step 3: Map reads to unitig sequences and get BAM files](#step-3-map-reads-to-unitig-sequences-and-get-bam-files)
  - [Step 4: Run CoverM to get coverage of unitig sequences](#step-4-run-coverm-to-get-coverage-of-unitig-sequences)
  - [Step 5: Scan unitig sequences for single-copy marker genes and PHROGs](#step-5-scan-unitig-sequences-for-single-copy-marker-genes-and-phrogs)
- [Phables Usage](#phables-usage)
  - [Phables options](#phables-options)
  - [Input](#input)
  - [Example usage](#example-usage)
  - [Output](#output)
  - [Snakemake pipeline](#snakemake-pipeline)
- [Reporting Issues](#reporting-issues)
  

## Setting up Phables

### Downloading Phables

You can clone the Phables repository to your machine.

```bash
git clone https://github.com/Vini2/phables.git
```

Now go into the `phables` folder using the command

```bash
cd phables/
```

### Installing Phables

We recommend that you use [`conda`](https://docs.conda.io/en/latest/) to install and run. Once you have installed `conda`, make sure you are in the `phables` folder. Now run the following commands to create a `conda` environment and activate it to run Phables.

```bash
conda env create -f environment.yml
conda activate phables
```

If you prefer to use `pip` instead of `conda`, you can run the following command to install Phables using `pip`. Make sure you are in the `phables` folder.

```bash
pip install .
```

### Setting up Gurobi

The MFD implementation uses the linear programming solver [Gurobi](https://www.gurobi.com/). The `phables` conda environment and pip setup already include Gurobi. To handle large models without any model size limitations, you have to activate the (academic) license and add the key using the following command.

```bash
grbgetkey <KEY>
```

Please refer to further instructions at [https://www.gurobi.com/academia/academic-program-and-licenses/](https://www.gurobi.com/academia/academic-program-and-licenses/).

### Test the setup

After setting up, run the following command to ensure that Phables is working.

```bash
phables --help
```

## Pre-processing

Phables requires the sequencing reads from viral metagenomic samples to be run through [Hecatomb](https://hecatomb.readthedocs.io/en/latest/) and be preprocessed. The following steps explain the preprocessing steps required to be carried out before running Phables.

Let's assume that the paired-end read files end as `<sample_name>_R1.fastq.gz` and `<sample_name>_R2.fastq.gz`. We can save the sample names to a file using the following command.

```bash
ls <reads_folder> | grep "R1.fastq.gz" | sed -e 's/_R1.fastq.gz//g'  > sample-names.txt
```

Now all the sample names will be saved to the file named `sample-names.txt`. We will use this file in the following steps.

### Step 1: Run read samples through Hecatomb

Phables requires the assembly output from [Hecatomb](https://hecatomb.readthedocs.io/en/latest/). First, run your sequencing read samples through Hecatomb using the following command.

```bash
hecatomb run --reads <reads_folder>
```

You can find the assembly output in `hecatomb.out/processing/assembly/CONTIG_DICTIONARY/FLYE/`. Phables will be using the `assembly_graph.gfa` and `assembly_info.txt` files.

### Step 2: Obtain unitig sequences from assembly graph

The `assembly_graph.gfa` file contains the unitig sequences which should be extracted. For this, you can use the `gfa2fasta.py` script found in `phables_utils/support/` as follows. This script will output the unitig sequences into a FASTA file named `edges.fasta`.

```bash
python gfa2fasta.py --graph assembly_graph.gfa --assembler flye --output <output_folder>
```

### Step 3: Map reads to unitig sequences and get BAM files

Use Minimap2 to map the reads of all the samples in `sample-names.txt` to unitigs in the `edges.fasta` file from step 2 and Samtools to index the BAM files.

```bash
mkdir bam_files

for f in `cat sample-names.txt`; do
  minimap2 -t 64 -N 5 -ax sr edges.fasta reads/"$f"_R1.fastq.gz reads/"$f"_R2.fastq.gz | samtools view -F 3584 -b --threads 64 > bam_files/"$f".bam
  samtools index bam_files/"$f".bam bam_files/"$f".bam.bai
```

### Step 4: Run CoverM to get coverage of unitig sequences

Use [CoverM](https://github.com/wwood/CoverM) to obtain the coverage of unitigs in the `edges.fasta` in each sample found in `sample-names.txt`.

```bash
mkdir coverage_rpkm

for f in `cat sample-names.txt`; do
  coverm contig -m rpkm -1 reads/"$f"_R1.fastq.gz -2 reads/"$f"_R2.fastq.gz -r edges.fasta -t 64 --output-file coverage_rpkm/"$f"_rpkm.tsv
```

This command will produce a coverage file for each sample. You can combine the coverage values of multiple samples into one file by running the `combine_cov.py` script found in `phables_utils/support/` as follows.

```bash
python combine_cov.py --covpath coverage_rpkm --output <output_path>
```

The output file will be `coverage.tsv` where each row represents a unitig and each column represents a sample.

### Step 5: Scan unitig sequences for single-copy marker genes and PHROGs

Scan unitigs for single-copy marker genes using [FragGeneScan](https://omics.informatics.indiana.edu/FragGeneScan/) and [HMMER](http://hmmer.org/). The  `marker.hmm` file for single-copy marker genes can be found [here](https://github.com/metagentools/MetaCoAG/tree/develop/metacoag_utils/auxiliary). Then scan the sequences for [PHROGs](https://phrogs.lmge.uca.fr). The PHROGs profile database for MMseqs can be downloaded from [here](https://phrogs.lmge.uca.fr/READMORE.php).

```bash
# SMG
run_FragGeneScan.pl -genome=edges.fasta -out=edges.fasta.frag -complete=0 -train=complete -thread=8 1>edges.fasta.frag.out 2>edges.fasta.frag.err
hmmsearch --domtblout edges.fasta.hmmout --cut_tc --cpu 8 marker.hmm edges.fasta.frag.faa 1>edges.fasta.hmmout.out 2> edges.fasta.hmmout.err

# PHROGs
mmseqs createdb edges.fasta target_seq
mmseqs search target_seq phrogs_profile_db results_mmseqs ./tmp -s 7
mmseqs createtsv target_seq phrogs_profile_db results_mmseqs phrog_annot.tsv
```

## Phables Usage

### Phables options

You can see the following command-line options of Phables using `python phables --help`.

```
Usage: phables [OPTIONS]

  Phables: Phage bubbles resolve bacteriophage genomes in viral metagenomic
  samples.

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
  -as, --alignscore FLOAT   minimum alignment score (%) for phrog annotations
  -si, --seqidentity FLOAT  minimum sequence identity for phrog annotations
  -o, --output PATH         path to the output folder  [required]
  --help                    Show this message and exit.
```

### Input

Phables takes the following files/paths as input.

* Assembly graph file (`assembly_graph.gfa`)
* Path information of contigs (`assembly_info.txt`)
* Coverage of unitig sequences in each sample in `.tsv` format (`coverage.tsv`) from preprocessing step 4
* HMMER result for scan of single-copy marker genes (`edges.fasta.hmmout`)
* MMseqs result for scan of PHROGs (`phrog_annot.tsv`)
* Path to BAM files (`bam_files/`)

### Example usage

```bash
phables -g assembly_graph.gfa -p assembly_info.txt -hm edges.fasta.hmmout -ph phrog_annot.tsv -c coverage.tsv -b bam_files/ -o /output/path/
```

### Output

The output from Phables will be as follows.

* `resolved_paths.fasta` containing the resolved genomes
* `resolved_phages` folder containing the resolved genomes in individual FASTA files
* `resolved_genome_info.txt` containing the path name, coverage, length, GC content and unitig order of the resolved genomes
* `resolved_edges.fasta` containing the unitigs that make up the resolved genomes
* `resolved_component_info.txt` containing the details of the phage bubbles resolved

### Snakemake pipeline
To Do

## Reporting Issues

Phables is still under testing. If you want to test (or break) Phables give it a try and report any issues and suggestions under [Phables Issues](https://github.com/Vini2/phables/issues).
