<p align="center">
  <img src="https://raw.githubusercontent.com/Vini2/phables/master/phables_logo.png" width="600" title="phables logo" alt="phables logo">
</p>

Phables: from fragmented assemblies to high-quality bacteriophage genomes
===============

[![DOI](https://img.shields.io/badge/Preprint_DOI-10.1101/2023.04.04.535632-blue)](https://doi.org/10.1101/2023.04.04.535632)
[![DOI](https://zenodo.org/badge/516191931.svg)](https://zenodo.org/badge/latestdoi/516191931)
![GitHub](https://img.shields.io/github/license/Vini2/phables)
[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![Conda](https://img.shields.io/conda/v/bioconda/phables)
![Conda](https://img.shields.io/conda/dn/bioconda/phables)
[![PyPI version](https://badge.fury.io/py/phables.svg)](https://badge.fury.io/py/phables)
[![Downloads](https://static.pepy.tech/badge/phables)](https://pepy.tech/project/phables)
[![CI](https://github.com/Vini2/phables/actions/workflows/testing.yml/badge.svg)](https://github.com/Vini2/phables/actions/workflows/testing.yml)
[![CodeQL](https://github.com/Vini2/phables/actions/workflows/codeql.yml/badge.svg)](https://github.com/Vini2/phables/actions/workflows/codeql.yml)
[![Documentation Status](https://readthedocs.org/projects/phables/badge/?version=latest)](https://phables.readthedocs.io/en/latest/?badge=latest)
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/Vini2/phables/develop?color=8a35da)

Phables is a tool developed to resolve bacteriophage genomes using phage bubbles in viral metagenomic data. It models cyclic phage-like components in the viral metagenomic assembly as flow networks, models as a minimum flow decomposition problem and resolves genomic paths corresponding to flow paths determined. Phables uses the [Minimum Flow Decomposition via Integer Linear Programming](https://github.com/algbio/MFD-ILP) implementation to obtain the flow paths.

For detailed instructions on installation and usage, please refer to the [**documentation hosted at Read the Docs**](https://phables.readthedocs.io/en/latest/).

**NEW:** Phables is now available on bioconda at [https://anaconda.org/bioconda/phables](https://anaconda.org/bioconda/phables) and on PyPI at [https://pypi.org/project/phables/](https://pypi.org/project/phables/). Feel free to pick your package manager, but we recommend that you use [`conda`](https://docs.conda.io/en/latest/).


## Setting up Phables

### Option 1: Installing Phables using conda (recommended)

You can install Phables from bioconda at [https://anaconda.org/bioconda/phables](https://anaconda.org/bioconda/phables). Make sure you have [`conda`](https://docs.conda.io/en/latest/) installed.

```bash
# create conda environment and install phables
conda create -n phables -c conda-forge -c anaconda -c bioconda phables

# activate environment
conda activate phables
```

Now you can go to [Setting up Gurobi](#setting-up-gurobi) to configure Gurobi.

### Option 2: Installing Phables using pip

You can install Phables from PyPI at [https://pypi.org/project/phables/](https://pypi.org/project/phables/). Make sure you have [`pip`](https://pip.pypa.io/en/stable/) and [`mamba`](https://mamba.readthedocs.io/en/latest/index.html) installed.

```bash
pip install phables
```

Now you can go to [Setting up Gurobi](#setting-up-gurobi) to configure Gurobi.

### Setting up Gurobi

The MFD implementation uses the linear programming solver [Gurobi](https://www.gurobi.com/). The `phables` conda environment and pip setup does not include Gurobi. You have to install Gurobi using one of the following commands depending on your package manager.

```bash
# conda
conda install -c gurobi gurobi

# pip
pip install gurobipy
```

To handle large models without any model size limitations, once you have installed Gurobi, you have to activate the (academic) license and add the key using the following command. You only have to do this once.

```bash
grbgetkey <KEY>
```

You can refer to further instructions at [https://www.gurobi.com/academia/academic-program-and-licenses/](https://www.gurobi.com/academia/academic-program-and-licenses/). 

### Test the installation

After setting up, run the following command to print out the Phables help message.

```bash
phables --help
```

## Quick Start Guide

Phables is powered by [Snaketool](https://github.com/beardymcjohnface/Snaketool) which packs in all the setup, testing, preprocessing and running steps into an easy-to-use pipeline.

### Setup the databases

```bash
# Download and setup the databases - you only have to do this once
phables install
```

### Run on test data

```bash
phables test
```

### Run on your own data

```bash
# Run Phables
# locally: using 8 threads (default is 1 thread)
phables run --input assembly_graph.gfa --reads fastq/ --threads 8
```

Please refer to the [**documentation hosted at Read the Docs**](https://phables.readthedocs.io/en/latest/) for further information on how to run Phables.


##  Issues and Questions

Phables is still under testing. If you want to test (or break) Phables give it a try and report any issues and suggestions under [Phables Issues](https://github.com/Vini2/phables/issues).

If you come across any questions, please have a look at the [Phables FAQ page](https://phables.readthedocs.io/en/latest/faq/). If your question is not here, feel free to post it under [Phables Issues](https://github.com/Vini2/phables/issues).


## Contributing to Phables

Are you interested in contributing to the Phables project? If so, you can check out the contributing guidelines in [CONTRIBUTING.md](https://github.com/Vini2/phables/blob/develop/CONTRIBUTING.md).


## Acknowledgement

Phables uses the [Gurobi](https://www.gurobi.com/) implementation of [MFD-ILP](https://github.com/algbio/MFD-ILP) and code snippets from [STRONG](https://github.com/chrisquince/STRONG), [METAMVGL](https://github.com/ZhangZhenmiao/METAMVGL), [GraphBin](https://github.com/metagentools/GraphBin), [MetaCoAG](https://github.com/metagentools/MetaCoAG) and [Hecatomb](https://hecatomb.readthedocs.io/en/latest/). Special thanks are owed to [Ryan Wick](https://github.com/rrwick) for developing [Bandage](https://rrwick.github.io/Bandage/) to visualise assembly graphs, which I heavily rely upon to investigate, develop and optimise my methods. The Phables logo was designed by [Amber Skye](https://fame.flinders.edu.au/people/2021/01/01/amber-cook).

## Citation
The Phables manuscript is currently under review and the preprint is available on biorxiv at DOI: [10.1101/2023.04.04.535632](https://www.biorxiv.org/content/10.1101/2023.04.04.535632v1). 

If you use Phables in your work, please cite Phables as,

> Vijini Mallawaarachchi, Michael J. Roach, Bhavya Papudeshi, Sarah K Giles, Susanna R Grigson, Przemyslaw Decewicz, George Bouras, Ryan D Hesse, Laura K Inglis, Abbey LK Hutton, Elizabeth A Dinsdale, Robert A Edwards. Phables: from fragmented assemblies to high-quality bacteriophage genomes. bioRxiv 2023.04.04.535632; doi: https://doi.org/10.1101/2023.04.04.535632.

```bibtex
@article{Mallawaarachchi2023.04.04.535632,
	author = {Vijini Mallawaarachchi and Michael J Roach and Bhavya Papudeshi and Sarah K Giles and Susanna R Grigson and Przemyslaw Decewicz and George Bouras and Ryan D Hesse and Laura K Inglis and Abbey LK Hutton and Elizabeth A Dinsdale and Robert A Edwards},
	title = {Phables: from fragmented assemblies to high-quality bacteriophage genomes},
	elocation-id = {2023.04.04.535632},
	year = {2023},
	doi = {10.1101/2023.04.04.535632},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Microbial communities found within the human gut have a strong influence on human health. Intestinal bacteria and viruses influence gastrointestinal diseases such as inflammatory bowel disease. Viruses infecting bacteria, known as bacteriophages, play a key role in modulating bacterial communities within the human gut. However, the identification and characterisation of novel bacteriophages remain a challenge. Available tools use similarities between sequences, nucleotide composition, and the presence of viral genes/proteins. Most available tools consider individual contigs to determine whether they are of viral origin. As a result of the challenges in viral assembly, fragmentation of viral genomes can occur, leading to the need for new approaches in viral identification. We introduce Phables, a new computational method to resolve bacteriophage genomes from fragmented viral metagenomic assemblies. Phables identifies bacteriophage-like components in the assembly graph, models each component as a flow network, and uses graph algorithms and flow decomposition techniques to identify genomic paths. Experimental results of viral metagenomic samples obtained from different environments show that over 80\% of the bacteriophage genomes resolved by Phables have high quality and are longer than the individual contigs identified by existing viral identification tools.},
	URL = {https://www.biorxiv.org/content/early/2023/04/04/2023.04.04.535632},
	eprint = {https://www.biorxiv.org/content/early/2023/04/04/2023.04.04.535632.full.pdf},
	journal = {bioRxiv}
}

```

Also, please cite the following tools/databases used by Phables.

* Roach MJ, Pierce-Ward NT, Suchecki R, Mallawaarachchi V, Papudeshi B, et al. Ten simple rules and a template for creating workflows-as-applications. PLOS Computational Biology 18(12) (2022): e1010705. [https://doi.org/10.1371/journal.pcbi.1010705](https://doi.org/10.1371/journal.pcbi.1010705)
* Terzian P, Olo Ndela E, Galiez C, Lossouarn J, Pérez Bucio RE, Mom R, Toussaint A, Petit MA, Enault F. PHROG: families of prokaryotic virus proteins clustered using remote homology. NAR Genomics and Bioinformatics, Volume 3, Issue 3, lqab067 (2021). [https://doi.org/10.1093/nargab/lqab067](https://doi.org/10.1093/nargab/lqab067)
* Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotechnol 35, 1026–1028 (2017). [https://doi.org/10.1038/nbt.3988](https://doi.org/10.1038/nbt.3988)
* Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100 (2018). [https://doi.org/10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191)
* Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools, Bioinformatics, Volume 25, Issue 16, Pages 2078–2079 (2009). [https://doi.org/10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352)
* Hagberg AA, Schult DA, and Swart PJ. Exploring network structure, dynamics, and function using NetworkX. In Proceedings of the 7th Python in Science Conference (SciPy2008), Gäel Varoquaux, Travis Vaught, and Jarrod Millman (Eds), (Pasadena, CA USA), pp. 11–15 (2008).
