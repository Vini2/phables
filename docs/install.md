# Setting up Phables

Phables is available on bioconda at [https://anaconda.org/bioconda/phables](https://anaconda.org/bioconda/phables) and on PyPI at [https://pypi.org/project/phables/](https://pypi.org/project/phables/). Feel free to pick your package manager, but we recommend that you use [`conda`](https://docs.conda.io/en/latest/).

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

## Setting up Gurobi

The MFD implementation uses the linear programming solver [Gurobi](https://www.gurobi.com/). We chose Gurobi over open source solvers as Gurobi is fast and can solve large models (check the performance comparison at [https://www.gurobi.com/resources/open-source-linear-and-mixed-integer-programming-software-and-solvers/](https://www.gurobi.com/resources/open-source-linear-and-mixed-integer-programming-software-and-solvers/)).

The `phables` conda environment and pip setup does not include Gurobi. You have to install Gurobi using one of the following commands depending on your package manager.

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

You can refer to further instructions at [https://www.gurobi.com/academia/academic-program-and-licenses/](https://www.gurobi.com/academia/academic-program-and-licenses/). Please note that this academic lisence is a file-based host-locked license, meaning that you can only run Gurobi on the machine that the license was obtained for. If you want to run on a cluster, you will have to contact your system admin and setup a floating network license. You can find more details at [https://support.gurobi.com/hc/en-us/articles/360013195412-How-do-I-obtain-a-free-academic-license-for-a-cluster-or-a-shared-computer-lab-](https://support.gurobi.com/hc/en-us/articles/360013195412-How-do-I-obtain-a-free-academic-license-for-a-cluster-or-a-shared-computer-lab-).

## Test the installation

After setting up, run `phables --help` to print out the Phables help message.

```bash
Usage: phables [OPTIONS] COMMAND [ARGS]...

  Phables: from fragmented assemblies to high-quality bacteriophage genomes.
  Please refer the full documentation available on Read the Docs at
  https://phables.readthedocs.io/

Options:
  -v, --version  Show the version and exit.
  -h, --help     Show this message and exit.

Commands:
  run       Run Phables
  install   Install databases
  test      Test Phables
  config    Copy the system default config file
  citation  Print the citation(s) for this tool
```

## Setup the databases

Now run the following command to download and setup the required databases.

```bash
phables install
```

## Run on the test data

Then run the following command to launch the test run and ensure that Phables is working.

```bash
phables test
```

If the test run completes without any issues, we are good to go.

## Build the docs

Optionally, the complete documentation of Phables including these pages can be built using [MkDocs](https://www.mkdocs.org/) as follows.

```bash
# install mkdocs
pip install mkdocs

# go to your installation directory
cd /path/to/phables

# build
mkdocs build
```