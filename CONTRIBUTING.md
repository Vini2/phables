# Contributing to Phables project

We love to have your contributions to the Phables project, whether it's:
* Reporting a bug
* Submitting a fix
* Proposing new features

## Clone and install Phables onto your machine

First, make sure you have [git](https://github.com/git-guides/install-git) installed on your machine.

On GitHub, [fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo) the Phables repository and [clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) it to your machine.

```bash
# clone repository to your local machine
git clone https://github.com/Vini2/phables.git
```

Move to the Phables directory 

```bash
cd phables
```

Create and activate the conda environment. Make sure to have [`conda`](https://docs.conda.io/en/latest/) installed.

```bash
# Create the phables environment
conda env create -f build/environment.yml

# Activate the phables environment
conda activate phables
```

Now install Phables via [`pip`](https://pip.pypa.io/en/stable/).

```bash
pip install -e .
```

## Test Phables installation

Print the help message using the following command.

```bash
phables -h
```

Use the following command to launch the phables test run and all the tests should pass.

```bash
phables test
```

## Coding Style

We adhere to the [PEP 8](https://peps.python.org/pep-0008/) style guide. 

Before committing, make sure to run [`black`](https://pypi.org/project/black/) and [`isort`](https://pypi.org/project/isort/).

## Report bugs using GitHub's issues

We use GitHub issues to track public bugs. Report a bug by opening a new issue in GitHub [issues](https://github.com/Vini2/phables/issues). You will get to select between templates for bug report and feature request. If none of these templates match what you want to report, you can use the custom issue template.

## Committing code

Once you have finished coding and all the tests pass, commit your code and make a pull request. 

```bash
# Add changed/added files
git add <file name>

# Commit changes
git commit -m "<commit message>"

# Push changes
git push
```

Make sure to follow the commit style of [c3dev](https://github.com/cogent3/c3dev/wiki#style-for-commit-messages). Relevant prefixes are replicated below for convenience.

| **Commit Prefix** | **For**                                       |
|-------------------|-----------------------------------------------|
| DEV:              | development tool or utility                   |
| DOC:              | documentation                                 |
| TST:              | addition or modification of tests             |
| REL:              | related to a release                          |
| MAINT:            | maintenance commit (refactoring, typos, etc.) |
| BUG:              | bug fix                                       |
| REV:              | revert an earlier commit                      |


Your contribution will be reviewed before accepting it. 

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

## References

This document was adapted from the open-source contribution guidelines for [Transcriptase](https://github.com/briandk/transcriptase-atom/blob/master/CONTRIBUTING.md) and [c3dev](https://github.com/cogent3/c3dev/wiki/How-to-Contribute-Code).