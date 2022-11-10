<p align="center">
  <img src="phables_logo.png" width="700" title="phables logo" alt="phables logo">
</p>

# Phables: Phage bubbles resolve bacteriophage genomes in viral metagenomic samples

[![CI](https://github.com/Vini2/phables/actions/workflows/testing.yml/badge.svg)](https://github.com/Vini2/phables/actions/workflows/testing.yml)
![GitHub](https://img.shields.io/github/license/Vini2/phables)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

`phables` is a tool developed to resolve bacteriophage genomes using phage bubbles in metagenomic data. It models cyclic phage-like components in the metagenomic assembly as flow networks and resolves genomic paths corresponding to flow paths determined. `phables` uses the [Minimum Flow Decomposition via  Integer Linear Programming](https://github.com/algbio/MFD-ILP) implementation to obtain the flow paths.

The figure below depicts the workflow of `phables`.

![](Phables_workflow.png)

`phables` requires the assembly output from [Hecatomb](https://hecatomb.readthedocs.io/en/latest/). 

## Setting up phables

### Downloading phables

You can clone the `phables` repository to your machine.

```
git clone https://github.com/Vini2/phables.git
```

Now go into the `phables` folder using the command

```
cd phables/
```

### Using `conda`

Once you have installed `conda`, make sure you are in the `phables` folder. Now run the following commands to create a `conda` environment and activate it to run `phables`.

```
conda env create -f environment.yml
conda activate phables
```

### Using `pip`
You can run the following command to install phables using `pip`. Make sure you are in the `phables` folder.

```
pip install .
```

### Setting up Gurobi

The MFD implementation uses the linear programming solver [Gurobi](https://www.gurobi.com/). The `phables` environment already comes with Gurobi installed. You have to activate the (academic) license and add the key using the following command.

```
grbgetkey <KEY>
```

Please refer to further instructions at [https://www.gurobi.com/academia/academic-program-and-licenses/](https://www.gurobi.com/academia/academic-program-and-licenses/).

### Test the setup

After setting up, run the following command to ensure that `phables` is working.

```
phables --help
```

## Usage

You can see the following command-line options of `phables` using `python phables --help`.

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

### Example usage

```
phables -g assembly_graph.gfa -p assembly_info.txt -hm edges.fasta.hmmout -ph phrog_annot.tsv -c coverage.tsv -b bam_files/ -o /output/path/
```

## Issues

`phables` is still under testing. Please report any issues and suggestions under [Phables Issues](https://github.com/Vini2/phables/issues).
