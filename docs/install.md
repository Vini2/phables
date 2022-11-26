# Setting up Phables

Phables is now available on PyPI at [https://pypi.org/project/phables/](https://pypi.org/project/phables/). If you directly installed Phables from [PyPI](https://pypi.org/project/phables/), you can skip the next two steps and go to [Setting up Gurobi](#setting-up-gurobi).

## Downloading Phables

You can clone the Phables repository to your machine.

```bash
git clone https://github.com/Vini2/phables.git
```

Now go into the `phables` folder using the command

```bash
cd phables/
```

## Installing Phables

We recommend that you use [`conda`](https://docs.conda.io/en/latest/) to install and run. Once you have installed `conda`, make sure you are in the `phables` folder. Now run the following commands to create a `conda` environment and activate it to run Phables.

```bash
conda env create -f environment.yml
conda activate phables
```

If you prefer to use `pip` instead of `conda`, you can run the following command to install Phables using `pip`. Make sure you are in the `phables` folder.

```bash
pip install .
```

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