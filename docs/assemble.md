# Assembly

Phables requires paired-end sequencing reads from viral metagenomic samples to be assembled. The following steps explain the steps required to be carried out beforehand.


## Paired-end reads files

Please make sure that the names of the read files are in the following format. Assuming that your paired-end sequencing reads are in the folder `fastq`, please make sure that the reads are in the format `{sampleName}_R1{fileExtension}`. `fileExtension` can be `.fq`, `.fastq`, `.fq.gz` or `.fastq.gz`.  For example,

```
sample1_R1.fastq.gz
sample1_R2.fastq.gz
sample2_R1.fastq.gz
sample2_R2.fastq.gz
...
```

## Assemble the samples

Phables requires the assembly graph file in **Graphical Fragment Assembly (GFA)** format. You can use any assembler that produces the assembly graph in GFA format to assemble your samples OR you can convert a FASTG file to GFA format.

If you have multiple samples, you can do a cross assembly or pool together reads and do a co-assembly.



Now we are ready to run Phables.