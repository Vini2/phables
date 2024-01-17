# Assembly

Phables requires paired-end sequencing reads from viral metagenomic samples to be assembled. The following steps explain the steps required to be carried out beforehand.


## Paired-end read files for short read assembly

Please make sure that the names of the read files are in the following format. Assuming that your paired-end sequencing reads are in the folder `fastq`, please make sure that the reads are in the format `{sampleName}{pattern}{fileExtension}`. `fileExtension` can be `.fq`, `.fastq`, `.fq.gz` or `.fastq.gz`.

Please make sure that your file `pattern` matches one of the following patterns.

```
_R1_ and _R2_
_R1. and _R2.
.R1. and .R2.
.R1_ and .R2_
_1_ and _2_
_1. and _2.
.1. and .2.
.1_ and .2_
```

For example, your read files can be

```
sample1_R1.fastq.gz
sample1_R2.fastq.gz
sample2_R1.fastq.gz
sample2_R2.fastq.gz
...
```

or

```
sample1_1.fq.gz
sample1_2.fq.gz
sample2_1.fq.gz
sample2_2.fq.gz
...
```

## Long read assemblies

If you are using long read datasets, there is no specific naming format for the read files.

## Assemble the samples

Phables requires the assembly graph file in **Graphical Fragment Assembly (GFA)** format. You can use any assembler that produces the assembly graph in GFA format to assemble your samples OR you can convert a FASTG file to GFA format.

If you have multiple samples you can pool together reads and do a co-assembly.

### Recommended assemblers and tools

* [MEGAHIT](https://github.com/voutcn/megahit) - for short reads
* [metaSPAdes](https://github.com/ablab/spades) - for short reads
* [metaFlye](https://github.com/fenderglass/Flye) - for long reads
* [Hecatomb](https://github.com/shandley/hecatomb) - viral analysis pipeline which produces a pooled assembly of the provided data

Now we are ready to run Phables.