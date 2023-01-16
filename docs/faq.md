# Frequently Asked Questions

## Q1: Where can I get help with issues?

If you come across any issues while using Phables, you can open an [issue on GitHub](https://github.com/Vini2/phables/issues) and we will look into it. Phables is still under development and testing, so we expect that there will still be bugs and unhandled exceptions in the code. 

## Q2: How can I run Phables on a cluster?

Phables can be run on a cluster as a normal package. However, your cluster should have a [floating network license](https://en.wikipedia.org/wiki/Floating_licensing) for Gurobi for the `run` subcommand to execute properly.

**Gurobi license for a cluster:** You will have to contact your system admin and setup a floating network license for the cluster. You can find more details at [https://support.gurobi.com/hc/en-us/articles/360013195412-How-do-I-obtain-a-free-academic-license-for-a-cluster-or-a-shared-computer-lab-](https://support.gurobi.com/hc/en-us/articles/360013195412-How-do-I-obtain-a-free-academic-license-for-a-cluster-or-a-shared-computer-lab-).

## Q3: Should I run Hecatomb before running Phables?

Yes, you have to run Hecatomb on your your raw reads and provide the `hecatomb.out` folder as inputs to the `phables run` command. Phables will scan through the `hecatomb.out` folder and automatically track down the assembly files.

## Q4: What can I do after running Phables?

Once you have run Phables, check out the [EVALUATION](https://phables.readthedocs.io/en/latest/quality/) section where you can read on how to check and compare the quality of the resolved genomes, interpret graph statistics and visualise the results.