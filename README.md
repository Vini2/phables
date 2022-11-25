<p align="center">
  <img src="https://raw.githubusercontent.com/Vini2/phables/master/phables_logo.png" width="600" title="phables logo" alt="phables logo">
</p>

Phables: Phage bubbles resolve bacteriophage genomes in viral metagenomic samples
===============

[![CI](https://github.com/Vini2/phables/actions/workflows/testing.yml/badge.svg)](https://github.com/Vini2/phables/actions/workflows/testing.yml)
![GitHub](https://img.shields.io/github/license/Vini2/phables)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Phables is a tool developed to resolve bacteriophage genomes using phage bubbles in viral metagenomic data. It models cyclic phage-like components in the viral metagenomic assembly as flow networks, models as a minimum flow decomposition problem and resolves genomic paths corresponding to flow paths determined. Phables uses the [Minimum Flow Decomposition via  Integer Linear Programming](https://github.com/algbio/MFD-ILP) implementation to obtain the flow paths.

For detailed instructions on installation and usage, please refer to the [**documentation hosted at Read the Docs**](https://phables.readthedocs.io/en/latest/).

**NEW:** Phables is now available on PyPI at [https://pypi.org/project/phables/](https://pypi.org/project/phables/).

<p align="center">
  <img src="https://raw.githubusercontent.com/Vini2/phables/master/Phables_workflow.png" title="phables workflow" alt="phables workflow">
</p>

## Setting up Phables

If you directly installed Phables from [PyPI](https://pypi.org/project/phables/), you can skip the next two steps and go to [Setting up Gurobi](#settingup-gurobi).

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

## Example usage

```bash
phables -g assembly_graph.gfa -p assembly_info.txt -hm edges.fasta.hmmout -ph phrog_annot.tsv -c coverage.tsv -b bam_files/ -o /output/path/
```

Please refer to the [**documentation hosted at Read the Docs**](https://phables.readthedocs.io/en/latest/) for further information on how to obtain/format the inputs.

## Reporting Issues

Phables is still under testing. If you want to test (or break) Phables give it a try and report any issues and suggestions under [Phables Issues](https://github.com/Vini2/phables/issues).


## Acknowledgement

Phables uses the Gurobi implementation of [MFD-ILP](https://github.com/algbio/MFD-ILP) and code snippets from [STRONG](https://github.com/chrisquince/STRONG), [METAMVGL](https://github.com/ZhangZhenmiao/METAMVGL), [GraphBin](https://github.com/metagentools/GraphBin) and [MetaCoAG](https://github.com/metagentools/MetaCoAG).