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

## Recommended assemblers and tools

### MEGAHIT

You can use [MEGAHIT](https://github.com/voutcn/megahit) to assemble your paired-end short read data.

```bash
megahit -1 reads_1.fastq -2 reads_2.fastq -o megahit_out
```

By default, MEGAHIT does not produce an assembly graph file. To get the assembly graph file, you have to to run `contig2fastg` command from the MEGAHIT toolkit to build the assembly graph. `contig2fastg` requires you to input the k-mer size used for the assembly. You can get the k-mer size from the contig IDs in the `final.contigs.fa` file. For example, you can use the `grep` command to print out the contig IDs as follows.

```bash
grep "^>" final.contigs.fa
```

For example, you will get an output as follows. Here the k-mer size is 141 as denoted by `k141`.

```bash
>k141_1456397 flag=0 multi=11.7570 len=1137
>k141_1235266 flag=0 multi=13.6963 len=1254
>k141_131192 flag=1 multi=47.8430 len=1510
>k141_1566081 flag=0 multi=9.6645 len=1372
...
```

Using the `k` value as 141, now you can run the `contig2fastg` command as follows.

```bash
megahit_toolkit contig2fastg 141 final.contigs.fa > final.graph.fastg
```

The MEGAHIT toolkit will result in a FASTG file which you can convert to GFA using [fastg2gfa](https://github.com/lh3/gfa1/blob/master/misc/fastg2gfa.c).

```bash
fastg2gfa final.graph.fastg > final.graph.gfa
```

If you want to run Phables on an assembly from a different `k` value output found in the MEGAHIT output folder `intermediate_contigs`, please make sure to build the `.fastg` file from the `.fa` file with the correct `k` value. For example, if you want to run Phables on the contigs from `k79.final.contigs.fa`, you should first build the corresponding `k79.final.graph.fastg` file and then run `fastg2gfa` as follows.

```bash
fastg2gfa 79.final.graph.fastg > 79.final.graph.gfa
```

### metaSPAdes

You can use [metaSPAdes](https://github.com/ablab/spades) to assemble your paired-end short read data. 

```bash
spades.py --meta -1 reads_1.fastq -2 reads_2.fastq -o metaspades_output -t 16
```

After the assembly finished, the output will contain the assembly graph file as `assembly_graph_after_simplification.gfa`.

### metaFlye

You can use [metaFlye](https://github.com/fenderglass/Flye) to assemble your long read data.

```bash
flye --meta --nano-raw reads.fasta --out-dir metaflye_output --threads 16
```

After the assembly finished, the output will contain the assembly graph file as `assembly_graph.gfa`.

### Hecatomb

You can use [Hecatomb](https://github.com/shandley/hecatomb) which is a viral analysis pipeline to obtain a pooled assembly of your short read or long read data contained in a folder named `reads`. You can run hecatomb as follows. Note that you only need to run the assembly module to process your data for Phables.

```bash
hecatomb run --reads reads/ assembly 
```

After the assembly finished, the output will contain the assembly graph file as `cross_assembly.gfa`.

Now we are ready to run Phables.