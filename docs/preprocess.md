# Pre-processing

Phables requires the sequencing reads from viral metagenomic samples to be run through [Hecatomb](https://hecatomb.readthedocs.io/en/latest/) and be preprocessed. The following steps explain the preprocessing steps required to be carried out before running Phables.

Let's assume that the paired-end read files end as `<sample_name>_R1.fastq.gz` and `<sample_name>_R2.fastq.gz`. We can save the sample names to a file using the following command.

```bash
ls <reads_folder> | grep "R1.fastq.gz" | sed -e 's/_R1.fastq.gz//g'  > sample-names.txt
```

Now all the sample names will be saved to the file named `sample-names.txt`. We will use this file in the following steps.

## Step 1: Run read samples through Hecatomb

Phables requires the assembly output from [Hecatomb](https://hecatomb.readthedocs.io/en/latest/). First, run your sequencing read samples through Hecatomb using the following command.

```bash
hecatomb run --reads <reads_folder>
```

You can find the assembly output in `hecatomb.out/processing/assembly/CONTIG_DICTIONARY/FLYE/`. Phables will be using the `assembly_graph.gfa` and `assembly_info.txt` files.

## Step 2: Obtain unitig sequences from assembly graph

The `assembly_graph.gfa` file contains the unitig sequences which should be extracted. For this, you can use the `gfa2fasta.py` script found in `phables_utils/support/` as follows. This script will output the unitig sequences into a FASTA file named `edges.fasta`.

```bash
python gfa2fasta.py --graph assembly_graph.gfa --assembler flye --output <output_folder>

# If you installed using pip
gfa2fasta --graph assembly_graph.gfa --assembler flye --output <output_folder>
```

## Step 3: Map reads to unitig sequences and get BAM files

Use Minimap2 to map the reads of all the samples in `sample-names.txt` to unitigs in the `edges.fasta` file from step 2 and Samtools to index the BAM files.

```bash
mkdir bam_files

for f in `cat sample-names.txt`
do
  minimap2 -t 64 -N 5 -ax sr edges.fasta reads/"$f"_R1.fastq.gz reads/"$f"_R2.fastq.gz | samtools view -F 3584 -b --threads 64 > bam_files/"$f".bam
  samtools index bam_files/"$f".bam bam_files/"$f".bam.bai
done
```

## Step 4: Run CoverM to get coverage of unitig sequences

Use [CoverM](https://github.com/wwood/CoverM) to obtain the coverage of unitigs in the `edges.fasta` in each sample found in `sample-names.txt`.

```bash
mkdir coverage_rpkm

for f in `cat sample-names.txt`
do
  coverm contig -m rpkm -1 reads/"$f"_R1.fastq.gz -2 reads/"$f"_R2.fastq.gz -r edges.fasta -t 64 --output-file coverage_rpkm/"$f"_rpkm.tsv
done
```

This command will produce a coverage file for each sample. You can combine the coverage values of multiple samples into one file by running the `combine_cov.py` script found in `phables_utils/support/` as follows.

```bash
python combine_cov.py --covpath coverage_rpkm --output <output_path>

# If you installed using pip
combine_cov --covpath coverage_rpkm --output <output_path>
```

The output file will be `coverage.tsv` where each row represents a unitig and each column represents a sample.

## Step 5: Scan unitig sequences for single-copy marker genes and PHROGs

Scan unitigs for single-copy marker genes using [FragGeneScan](https://omics.informatics.indiana.edu/FragGeneScan/) and [HMMER](http://hmmer.org/). The  `marker.hmm` file for single-copy marker genes can be found [here](https://github.com/metagentools/MetaCoAG/tree/develop/metacoag_utils/auxiliary). Then scan the sequences for [PHROGs](https://phrogs.lmge.uca.fr). The PHROGs profile database for MMseqs can be downloaded from [here](https://phrogs.lmge.uca.fr/READMORE.php).

```bash
# SMG
run_FragGeneScan.pl -genome=edges.fasta -out=edges.fasta.frag -complete=0 -train=complete -thread=8 1>edges.fasta.frag.out 2>edges.fasta.frag.err
hmmsearch --domtblout edges.fasta.hmmout --cut_tc --cpu 8 marker.hmm edges.fasta.frag.faa 1>edges.fasta.hmmout.out 2> edges.fasta.hmmout.err

# PHROGs
mmseqs createdb edges.fasta target_seq
mmseqs search target_seq phrogs_profile_db results_mmseqs ./tmp -s 7
mmseqs createtsv target_seq phrogs_profile_db results_mmseqs phrog_annot.tsv --full-header
```

**Note**: If you face problems while running your sequences against the PHROGs database using the latest MMseqs version, try running with MMseqs version 13 or lower.