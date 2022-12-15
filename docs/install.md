# Setting up Phables

Phables is available on Anaconda.org at [https://anaconda.org/vijinim/phables](https://anaconda.org/vijinim/phables) and on PyPI at [https://pypi.org/project/phables/](https://pypi.org/project/phables/). Feel free to pick your package manager.

### Option 1: Installing Phables from Anaconda.org

You can install Phables from Anaconda.org at [https://anaconda.org/vijinim/phables](https://anaconda.org/vijinim/phables).

```bash
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

```bash
pip install phables
```

Now you can go to [Setting up Gurobi](#setting-up-gurobi) to configure Gurobi.

## Setting up Gurobi

The MFD implementation uses the linear programming solver [Gurobi](https://www.gurobi.com/). We chose Gurobi over open source solvers as Gurobi is fast and can solve large models (check the performance comparison at [https://www.gurobi.com/resources/open-source-linear-and-mixed-integer-programming-software-and-solvers/](https://www.gurobi.com/resources/open-source-linear-and-mixed-integer-programming-software-and-solvers/)).

The `phables` conda environment and pip setup already include Gurobi. To handle large models without any model size limitations, you have to activate the (academic) license and add the key using the following command.

```bash
grbgetkey <KEY>
```

You can refer to further instructions at [https://www.gurobi.com/academia/academic-program-and-licenses/](https://www.gurobi.com/academia/academic-program-and-licenses/). Please note that this academic lisence is a file-based host-locked license, meaning that you can only run Gurobi on the machine that the license was obtained for. If you want to run on a cluster, you will have to contact your system admin and setup a floating network license. You can find more details at [https://support.gurobi.com/hc/en-us/articles/360013195412-How-do-I-obtain-a-free-academic-license-for-a-cluster-or-a-shared-computer-lab-](https://support.gurobi.com/hc/en-us/articles/360013195412-How-do-I-obtain-a-free-academic-license-for-a-cluster-or-a-shared-computer-lab-).

## Test the setup

After setting up, run the following command to ensure that Phables is working.

```bash
phables --help
```
