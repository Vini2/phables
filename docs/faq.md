# Frequently Asked Questions

## General FAQs

### Q1: Where can I get help with issues?

If you come across any issues while using Phables, you can open an [issue on GitHub](https://github.com/Vini2/phables/issues) and we will look into it. Phables is still under development and testing, so we expect that there will still be bugs and unhandled exceptions in the code. 

### Q2: Can I use the assembly graph from any assembler?

Phables supports any assembly graph in GFA (`.gfa`) format. You can use any assembler that produces the assembly graph in GFA format to assemble your samples OR you can convert an assembly graph in FASTG format to GFA format using a tool such as [fastg2gfa](https://github.com/lh3/gfa1/blob/master/misc/fastg2gfa.c).

If you use metaSPAdes for assembly, you can use the `assembly_graph_after_simplification.gfa` as input for Phables.

### Q4: What can I do after running Phables?

Once you have run Phables, check out the [EVALUATION](https://phables.readthedocs.io/en/latest/quality/) section where you can read on how to check and compare the quality of the resolved genomes, interpret graph statistics and visualise the results.

### Q5: How can I find out which contigs were included in the resolved phages?

The `resolved_genome_info.txt` file contains the order of sequences in the assembly graph that were used to construct the genomes (refer to the example below). 

| Path                 | Case  | Coverage | Length | GC content        | Node order                                                                   |
|----------------------|-------|----------|--------|-------------------|------------------------------------------------------------------------------|
| phage_comp_0_cycle_1 | case3 | 644      | 45659  | 34.86059703453864 | ['49-', '5524+', '24979-', '5556+', '55+', '5540-', '67+', '4490+', '5554-'] |
| phage_comp_0_cycle_2 | case3 | 625      | 43427  | 35.03810993160937 | ['49-', '5522+', '24979-', '5558+', '55+', '65+', '67+', '4490+', '5498-']   |
| ...                  |       |          |        |                   |                                                                              |

The `Node order` column denotes the segment IDs from the assembly graph.

The mapping between the contigs and assembly graph segments depends on the assembler you use. 

* If you use MEGAHIT, the segments in the assembly graph are the contigs themselves. You can directly relate the `Node order` information as the contigs that make the paths.
* If you use an assembler such as SPAdes or Flye, the sequences represented in the assembly graph are **unitigs**, which make up contigs. The information on which unitigs make up the contigs can be found in, for example, `contigs.paths` file in SPAdes and `assembly_info.txt` file in Flye.

### Q6: Can I run Phables on mixed-microbial communities?

Phables was originally designed to run on viromic data, but it can also be used to study mixed-microbial communities. However, the current implementation of Phables filters any component with at least a single unitig encoding any bacterial single-copy marker gene and hence, prophages might be omitted in the final result. Also, some plasmids or [phage-plasmids](https://doi.org/10.1128/mbio.01851-22), can be identified by Phables as phages. Hence, users should perform further downstream analysis to ensure that the predicted genomes are actual phages. One option is to use a tool such as [PPR-Meta](https://github.com/zhenchengfang/PPR-Meta) to classify the genomes resolved from Phables into phages and plasmids.

### Q7: Can Phables identify prophages?

If a prophage is active, excises from the genome, and is replicating, Phables would identify it. However, if it is a cryptic prophage, Phables would not identify it as it will be integrated into the host genome and can be part of a larger bacterial component in the assembly graph. As Phables discards components having bacterial single-copy marker genes, such prophages will not be identified. 

Users can use specific tools to either identify prophages in bacterial genomes such as [Phispy](https://academic.oup.com/nar/article/40/16/e126/1027055) or [hafeZ](https://www.biorxiv.org/content/10.1101/2021.07.21.453177v1) or validate recovered prophage sequences from host-genomes in metagenomic sequences such as [CheckV](https://www.nature.com/articles/s41587-020-00774-7).


## Gurobi FAQs

### Q1: Gurobi installation conflicts and `grbgetkey` fails to run

If you come across conflicts when installing Gurobi in the `phables` environment and could not run the `grbgetkey` command properly, please follow the steps given below.

```bash
# Deactivate the phables environment
conda deactivate

# Remove phables environemnt
conda remove -n phables --all

# Create conda environment with phables and gurobi
conda create -n phables -c conda-forge -c anaconda -c bioconda -c gurobi phables gurobi
```

### Q2: Model too large for size-limited license

If you get the following error when running Phables, this means that you don't have a proper license to handle large models. 

```bash
Error code 10010: Model too large for size-limited license; visit https://www.gurobi.com/free-trial for a full license
```

You should get an academic license which is provided free of charge to your institutional email address. You can refer to further instructions at [https://www.gurobi.com/academia/academic-program-and-licenses/](https://www.gurobi.com/academia/academic-program-and-licenses/).

### Q3: HostID mismatch

If you get the following error when running Phables as a job on a cluster, you cannot use your academic license which is a file-based host-locked license, meaning that you can only run Gurobi on the machine that the license was obtained for.

```bash
Failed to set up a license
Error 10009: HostID mismatch (licensed to <host_1>, hostid is <host_2>)
```

You will have to contact your system admin and setup a floating network license. You can find more details at [https://support.gurobi.com/hc/en-us/articles/360013195412-How-do-I-obtain-a-free-academic-license-for-a-cluster-or-a-shared-computer-lab-](https://support.gurobi.com/hc/en-us/articles/360013195412-How-do-I-obtain-a-free-academic-license-for-a-cluster-or-a-shared-computer-lab-).

### Q4: License not valid for Gurobi version x

If you get the following error when running Phables, this means that the version in your license does not match the installed version. You can install the correct version of Gurobi to match your license or you can get a new license for the latest version installed.

```bash
ERROR - Error code 10009: Request denied: license not valid for Gurobi version 11
```

### Q5: How can I get a Gurobi license for a cluster?

If you want to run Phables on a cluster, your cluster should have a [floating network license](https://en.wikipedia.org/wiki/Floating_licensing) for Gurobi for the `run` subcommand to execute properly.

**Gurobi license for a cluster:** You will have to contact your system admin and setup a floating network license for the cluster. You can find more details at [https://support.gurobi.com/hc/en-us/articles/360013195412-How-do-I-obtain-a-free-academic-license-for-a-cluster-or-a-shared-computer-lab-](https://support.gurobi.com/hc/en-us/articles/360013195412-How-do-I-obtain-a-free-academic-license-for-a-cluster-or-a-shared-computer-lab-).

If your cluster has Gurobi already installed with the license setup, you can load the module as follows, prior to running Phables.

```bash
module load gurobi
```