![](images/phables_logo.png)

# Phage bubbles resolve bacteriophage genomes in viral metagenomic samples

Phables is a tool developed to resolve bacteriophage genomes using phage bubbles in viral metagenomic data. 
It models cyclic phage-like components in a viral metagenomic assembly graph as flow networks, models as a 
minimum flow decomposition problem and resolves genomic paths corresponding to flow paths determined. 
Phables uses the [Minimum Flow Decomposition via  Integer Linear 
Programming](https://github.com/algbio/MFD-ILP) implementation to obtain the flow paths.

## Motivation

Existing viral identification tools run contigs through a pre-trained model and predict whether or not they are of viral origin. However, contigs do not necessarily represent complete genomes as viral assemblies are not always perfect. Most of the existing metagenomic binning tools are optimised for bacterial metagenomes and cannot handle viral metagenomes efficiently.

We observed circular components in viral metagenome assembly graphs as shown below (visualisations obtained from [Bandage](https://rrwick.github.io/Bandage/)), suggesting that viral genomes are fragmented and variants exist.

![](images/components.png)

Phables was developed to recover phage-like components called "phage bubbles" that represent one or more bacteriophage genomes and resolve phage bubbles to obtain complete and high-quality genomes.

## Workflow

![](images/Phables_workflow.png)

