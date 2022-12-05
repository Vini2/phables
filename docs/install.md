# Setting up Phables

Phables is available on Anaconda.org at [https://anaconda.org/vijinim/phables](https://anaconda.org/vijinim/phables) and on PyPI at [https://pypi.org/project/phables/](https://pypi.org/project/phables/). Feel free to pick your package manager.

### Option 1: Installing Phables from Anaconda.org

You can install Phables from Anaconda.org at [https://anaconda.org/vijinim/phables](https://anaconda.org/vijinim/phables).

```
# add channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels gurobi

# create conda environment and install phables
conda create -n phables -c vijinim phables

# activate environment
conda activate phables
```

Now you can go to [Setting up Gurobi](#setting-up-gurobi) to configure Gurobi.

### Option 2: Installing Phables from PyPi

You can install Phables from PyPI at [https://pypi.org/project/phables/](https://pypi.org/project/phables/).

```
pip install phables
```

Now you can go to [Setting up Gurobi](#setting-up-gurobi) to configure Gurobi.

## Setting up Gurobi

The MFD implementation uses the linear programming solver [Gurobi](https://www.gurobi.com/). The `phables` conda environment and pip setup already include Gurobi. To handle large models without any model size limitations, you have to activate the (academic) license and add the key using the following command.

```bash
grbgetkey <KEY>
```

Please refer to further instructions at [https://www.gurobi.com/academia/academic-program-and-licenses/](https://www.gurobi.com/academia/academic-program-and-licenses/).

## Test the setup

After setting up, run the following command to ensure that Phables is working.

```bash
phables --help
```