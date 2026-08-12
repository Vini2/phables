# Phables Usage

Phables run options can be found using the `phables run -h` command.

```
Usage: phables run [OPTIONS] [SNAKE_ARGS]...

  Run Phables

Options:
  --output PATH                   Output directory  [default: phables.out]
  --configfile TEXT               Custom config file [default:
                                  (outputDir)/config.yaml]
  --threads INTEGER               Number of threads to use  [default: 1]
  --databases PATH                Path to databases directory [default:
                                  <install_dir>/databases, i.e. wherever
                                  `phables install` put them]
  --use-conda / --no-use-conda    Use conda for Snakemake rules  [default:
                                  use-conda]
  --conda-prefix PATH             Custom conda env directory
  --container PATH                container image with every per-rule tool
                                  already installed (container/Dockerfile),
                                  replacing --use-conda's per-rule env
                                  creation for the WHOLE workflow -- not just
                                  predict_3di (--prostt5-container). Needs
                                  --use-singularity passed as a trailing
                                  snakemake arg; don't also pass --use-conda
                                  alongside this.
  --profile TEXT                  Snakemake profile
  --snake-default TEXT            Customise Snakemake runtime args  [default:
                                  --rerun-incomplete, --printshellcmds,
                                  --nolock, --show-failed-logs]
  --input PATH                    Path to assembly graph file in .GFA format
                                  [required]
  --reads PATH                    Path to directory containing paired-end
                                  reads  [required]
  --minlength INTEGER             minimum length of circular unitigs to
                                  consider as standalone genomes, and the
                                  minimum length for any resolved LINEAR path
                                  to be kept -- circular paths are exempt
                                  regardless of length, since a closed cycle
                                  is itself strong completeness evidence
                                  [default: 2000]
  --mincov INTEGER                minimum coverage of paths to output
                                  [default: 10]
  --compcount INTEGER             maximum unitig count to consider a component
                                  [default: 200]
  --maxpaths INTEGER              maximum number of paths to resolve for a
                                  component  [default: 10]
  --mgfrac FLOAT                  length threshold to consider single copy
                                  marker genes  [default: 0.2]
  --build-tree / --no-build-tree  align resolved genomes (MAFFT) and build a
                                  phylogenetic tree (IQ-TREE via piqtree). Off
                                  by default -- slow at metagenome scale
                                  (thousands of genomes); run phylogenetics as
                                  a separate step on the resolved genomes
                                  instead unless you specifically need a per-
                                  run tree  [default: no-build-tree]
  --genecaller [pyrodigal-gv|fraggenescan]
                                  gene caller to use for unitig gene
                                  prediction  [default: pyrodigal-gv]
  --phagedetection [mmseqs|prostt5-foldseek]
                                  phage-gene detection method: mmseqs (PHROGs,
                                  default) or prostt5-foldseek (structural;
                                  needs --hallmark-db, --hallmark-categories
                                  and --prostt5-checkpoint)  [default: mmseqs]
  --gpu-backend [cpu|cuda|rocm]   PyTorch build for the ProstT5 conda env:
                                  cpu, cuda, or rocm (e.g. Setonix's MI250X
                                  nodes). Independent of --prostt5-cpu, which
                                  forces ProstT5 onto the CPU device at
                                  runtime even inside a GPU-capable env --
                                  this controls which env gets built
                                  [default: cpu]
  --foldseek-gpu                  use foldseek's CUDA GPU search mode for the
                                  hallmark scan (requires --gpu-backend cuda,
                                  a CUDA-capable foldseek build on PATH -- not
                                  the plain bioconda package -- and a *_gpu-
                                  suffixed, makepaddedseqdb-prepared hallmark
                                  DB)
  --hallmark-db PATH              foldseek hallmark structure subDB prefix,
                                  only used with --phagedetection
                                  prostt5-foldseek. Auto-resolves to the copy
                                  `phables install` fetches under --databases
                                  (docs/hallmark_db.md) when not set --
                                  override only if you built your own with
                                  build_hallmark_db.py
  --hallmark-categories PATH      hallmark PHROG categories TSV, only used
                                  with --phagedetection prostt5-foldseek.
                                  Auto-resolves alongside --hallmark-db when
                                  not set -- see
                                  build_hallmark_db.py/docs/hallmark_db.md
  --hallmark-evalue FLOAT         maximum e-value for hallmark structural hits
                                  [default: 1e-08]
  --hallmark-minbits INTEGER      minimum bitscore for hallmark structural
                                  hits  [default: 0]
  --prostt5-checkpoint PATH       ProstT5 CNN prediction-head checkpoint
                                  (required for --phagedetection
                                  prostt5-foldseek)
  --prostt5-model TEXT            ProstT5 HuggingFace model identifier
                                  [default: Rostlab/ProstT5_fp16]
  --prostt5-model-dir PATH        directory to cache the ProstT5 model in
                                  [default: ~/.cache/prostt5]
  --prostt5-half-precision / --prostt5-full-precision
                                  run ProstT5 in half precision (ignored on
                                  CPU)  [default: prostt5-half-precision]
  --prostt5-cpu                   force ProstT5 onto CPU even if a GPU is
                                  available
  --prostt5-max-residues INTEGER  max total residues per ProstT5 batch --
                                  device-specific, tune per GPU  [default:
                                  4000]
  --prostt5-max-seq-len INTEGER   sequences longer than this flush a ProstT5
                                  batch immediately  [default: 4000]
  --prostt5-max-batch INTEGER     max sequences per ProstT5 batch -- device-
                                  specific, tune per GPU  [default: 20]
  --prostt5-container TEXT        container image with pholdlib + torch
                                  already installed (e.g. phold's own image),
                                  used instead of a conda env for predict_3di.
                                  Needs --use-singularity passed as a trailing
                                  snakemake arg -- --use-conda alone won't
                                  honour it. Overrides --gpu-backend for this
                                  rule.
  --evalue FLOAT                  maximum e-value for phrog annotations
                                  [default: 1e-10]
  --seqidentity FLOAT             minimum sequence identity for phrog
                                  annotations  [default: 0.3]
  --covtol INTEGER                coverage tolerance for extending subpaths
                                  [default: 100]
  --alpha FLOAT                   coverage multiplier for flow interval
                                  modelling  [default: 1.2]
  --longreads                     provide long reads as input (else defaults
                                  to short reads)
  --prefix TEXT                   prefix for genome identifier
  -h, --help                      Show this message and exit.

  
  If you use Phables in your work, please cite Phables as,
  
  Vijini Mallawaarachchi, Michael J Roach, Przemyslaw Decewicz, 
  Bhavya Papudeshi, Sarah K Giles, Susanna R Grigson, George Bouras, 
  Ryan D Hesse, Laura K Inglis, Abbey L K Hutton, Elizabeth A Dinsdale, 
  Robert A Edwards, Phables: from fragmented assemblies to high-quality 
  bacteriophage genomes, Bioinformatics, Volume 39, Issue 10, 
  October 2023, btad586, https://doi.org/10.1093/bioinformatics/btad586
  
  
  For more information on Phables please visit:
  https://phables.readthedocs.io/
  
  
  CLUSTER EXECUTION:
  phables run ... --profile [profile]
  For information on Snakemake profiles see:
  https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
  
  RUN EXAMPLES:
  Required:           phables run --input [assembly graph file]
  Specify threads:    phables run ... --threads [threads]
  Disable conda:      phables run ... --no-use-conda 
  Change defaults:    phables run ... --snake-default="-k --nolock"
  Add Snakemake args: phables run ... --dry-run --keep-going --touch
  Specify targets:    phables run ... print_stages
  Available targets:
      all             Run everything (default)
      preprocess      Run preprocessing only
      phables         Run phables (and preprocessing if needed)
      postprocess     Run postprocessing (with preprocessing and phables if needed)
      print_stages    List available stages
```

## Run options explained

* `--input` - assembly graph file in .GFA format
* `--reads` - folder containing paired-end read files
* `--minlength` - minimum length of circular unitigs to consider as standalone genomes, and the minimum length for any resolved *linear* path to be kept -- circular paths are exempt regardless of length, since a closed cycle is itself strong completeness evidence [default: 2000]
* `--mincov` - minimum coverage of paths to output [default: 10]
* `--compcount` - maximum unitig count to consider a component [default: 200]
* `--maxpaths` - maximum number of paths to resolve for a component [default: 10]
* `--mgfrac` - length threshold to consider single copy marker genes [default: 0.2]
* `--build-tree` / `--no-build-tree` - align resolved genomes (MAFFT) and build a phylogenetic tree (IQ-TREE via piqtree) as part of the run. Off by default -- slow at metagenome scale (thousands of genomes); run phylogenetics as a separate step on `resolved_paths.fasta` instead unless you specifically need a per-run tree [default: `no-build-tree`]
* `--genecaller` - gene caller for unitig gene prediction: `pyrodigal-gv` or `fraggenescan` [default: `pyrodigal-gv`]
* `--evalue` - maximum e-value for phrog annotations [default: 1e-10]
* `--seqidentity` - minimum sequence identity for phrog annotations [default: 0.3]
* `--covtol` - coverage tolerance for extending subpaths [default: 100]
* `--alpha` - coverage multiplier for flow interval modelling [default: 1.2]
* `--longreads` - provide long reads as input. If this flag is not provided phables defaults to short reads
* `--prefix` - prefix for genome identifier [default: None]
* `--output` - path to the output directory [default: `phables.out`]
* `--configfile` - custom config file [default: `(outputDir)/config.yaml`]
* `--threads` - number of threads to use  [default: 1]
* `--databases` - path to the databases directory [default: wherever `phables install` put them]
* `--use-conda` / `--no-use-conda` - use conda for Snakemake rules  [default: `use-conda`]
* `--conda-prefix` - custom conda env directory
* `--container` - run the whole workflow from a single container image instead of per-rule conda envs (needs `--use-singularity`, and `--no-use-conda`) -- see [Running from a single container](container.md)
* `--snake-default` - customise Snakemake runtime args  [default: `--rerun-incomplete, --printshellcmds, --nolock, --show-failed-logs`]

### Phage-gene detection: `--phagedetection`

Two independent methods for finding phage-like genes on unitigs, feeding the same downstream component-filtering/MFD logic either way:

* `mmseqs` (default) - PHROGs profile search via MMseqs2 against the database `phables install` fetches. No extra flags needed.
* `prostt5-foldseek` - structural detection instead of sequence homology: predicts each protein's 3Di structural alphabet with [ProstT5](https://github.com/mheinzinger/ProstT5), then searches it with Foldseek against a subDB of phage hallmark protein structures (see [Building the hallmark database](hallmark_db.md)). Useful for distant homologs sequence search misses, at the cost of needing a GPU to run at any reasonable speed (CPU works, just slowly).

`--phagedetection prostt5-foldseek` options:

* `--hallmark-db` / `--hallmark-categories` - auto-resolve to the copy `phables install` fetches under `--databases`; only pass these if you built your own with `build_hallmark_db.py`
* `--hallmark-evalue` - maximum e-value for hallmark structural hits [default: 1e-08]
* `--hallmark-minbits` - minimum bitscore for hallmark structural hits [default: 0]
* `--prostt5-checkpoint` - **required**: the ProstT5 CNN prediction-head checkpoint (phold's own `model.pt`)
* `--prostt5-model` - ProstT5 HuggingFace model identifier [default: `Rostlab/ProstT5_fp16`]
* `--prostt5-model-dir` - directory to cache the downloaded ProstT5 model in [default: `~/.cache/prostt5`]
* `--prostt5-half-precision` / `--prostt5-full-precision` - run ProstT5 in half precision (ignored on CPU) [default: `prostt5-half-precision`]
* `--prostt5-cpu` - force ProstT5 onto CPU even when a GPU is available
* `--prostt5-max-residues`, `--prostt5-max-seq-len`, `--prostt5-max-batch` - ProstT5 batching knobs, device-specific -- tune per GPU rather than trusting the defaults on unfamiliar hardware [defaults: 4000, 4000, 20]
* `--gpu-backend` - which PyTorch build the ProstT5 conda env solves against: `cpu`, `cuda`, or `rocm` (e.g. Setonix's MI250X). Independent of `--prostt5-cpu`, which forces the CPU device at runtime even inside a GPU-capable env -- this controls which env gets *built* [default: `cpu`]
* `--foldseek-gpu` - use Foldseek's CUDA GPU search mode for the hallmark scan. Requires `--gpu-backend cuda`, a CUDA-capable Foldseek build on `PATH` (not the plain bioconda package), and a `*_gpu`-suffixed, `makepaddedseqdb`-prepared hallmark DB
* `--prostt5-container` - use a container with `pholdlib`+torch already installed (e.g. phold's own image) instead of a conda env for the ProstT5 step. Needs `--use-singularity` passed as a trailing Snakemake arg -- `--use-conda` alone won't honour it. Overrides `--gpu-backend` for this rule.

## Example usage

Assuming your assembly graph file is `assembly_graph.gfa` and reads folder as `fastq`, you can run `phables` as follows.

### Using short reads

```bash
# Preprocess data using 8 threads (default is 1 thread)
phables run --input assembly_graph.gfa --reads fastq --threads 8
```

### Using long reads

```bash
# Preprocess data using 8 threads (default is 1 thread)
phables run --input assembly_graph.gfa --reads fastq --threads 8 --longreads
```

### Using structural (ProstT5+Foldseek) phage-gene detection

Same as above, plus `--phagedetection prostt5-foldseek` and its required checkpoint
(`--hallmark-db`/`--hallmark-categories` auto-resolve once `phables install` has
fetched them -- see [Building the hallmark database](hallmark_db.md)):

```bash
phables run --input assembly_graph.gfa --reads fastq --threads 8 \
    --phagedetection prostt5-foldseek \
    --prostt5-checkpoint /path/to/model.pt
```

On a GPU node, add `--gpu-backend cuda` (or `rocm`) so the ProstT5 conda env
actually builds against a GPU-capable torch -- without it, ProstT5 still runs,
just on CPU, which is slow for anything beyond a handful of proteins:

```bash
phables run --input assembly_graph.gfa --reads fastq --threads 8 \
    --phagedetection prostt5-foldseek \
    --prostt5-checkpoint /path/to/model.pt \
    --gpu-backend cuda
```

Note that you should provide the path to the GFA file to the `--input` parameter and the folder containing your sequencing reads to the `--reads` parameter. 

The output of Phables is set by default to `phables.out`. You can update the output path using the `--output` parameter for `phables run` as follows.

```bash
# Preprocess data using 8 threads (default is 1 thread)
phables run --input assembly_graph.gfa --reads fastq --output my_output_folder --threads 8
```

The `phables run` command will run preprocessing steps, perform genome resolution and the perform postprocessing steps.

## Output

Following is the folder structure of the Phables complete run.

```
phable.out
├── config.yaml  # config file
├── logs         # all log files
├── phables      # final phables results
├── phables.log  # phables master log
├── postprocess  # postprocessing results
└── preprocess   # preprocessing results
```

Phables will create 3 main folders `preprocess`, `phables` and `postprocess` for the different stages of execution.

### 1. `preprocess` - preprocessing results

The following preprocessing steps will be run and their corresponding files and folders can be found in the `preprocess` folder.

* Obtain unitig sequences from assembly graph - `edges.fasta`
* Call genes on unitigs (`--genecaller pyrodigal-gv` or `fraggenescan`) - `edges.fasta.frag.faa`
* Map reads to unitig sequences and get BAM files (via [CoverM](https://github.com/wwood/CoverM)) - `temp/*.bam` and `temp/*.bai`
* Calculate coverage of unitig sequences - `coverage.tsv`
* Scan unitig sequences for single-copy marker genes - `edges.fasta.hmmout`
* Scan unitig sequences for phage-like genes, via one of two methods depending on `--phagedetection`:
    * `mmseqs` (default): search against Prokaryotic Virus Remote Homologous Groups ([PHROGs](https://phrogs.lmge.uca.fr/)) - `phrogs_annotations.tsv`
    * `prostt5-foldseek`: ProstT5-predicted 3Di structures searched against the hallmark structure subDB (see [Building the hallmark database](hallmark_db.md)) - `hallmark/proteins_3di.fasta` and `hallmark_hits.tsv`

### 2. `phables` - genome resolution results

The following files and folders can be found inside the `phables` folder which are the main outputs of Phables.

* `resolved_paths.fasta` containing the resolved genomes
* `resolved_phages` folder containing the resolved genomes in individual FASTA files
* `resolved_genome_info.txt` containing the path name, coverage, length, GC content and unitig order of the resolved genomes
* `resolved_edges.fasta` containing the unitigs that make up the resolved genomes
* `unresolved_phage_like_edges.fasta` containing all the unresolved phage-like unitigs
* `all_phage_like_edges.fasta` containing sequences from all the phage-like components (both resolved and unresolved)
* `resolved_component_info.txt` containing the details of the phage bubbles resolved
* `component_phrogs.txt` containing PHROGs found in each component

### 3. `postprocess` - postprocessing results

The following postprocessing steps will be run and their corresponding files and folders can be found in the `postprocess` folder.

* Combine resolved genomes and unresolved edges - `genomes_and_unresolved_edges.fasta`
* Obtain read counts for resolved genomes and unresolved edges - `sample_genome_read_counts.tsv`
* Obtain mean coverage of resolved genomes and unresolved edges - `sample_genome_mean_coverage.tsv`
* Obtain RPKM coverage of resolved genomes and unresolved edges - `sample_genome_rpkm.tsv`
* (optional, `--build-tree`) Generate a multiple sequence alignment of the resolved genomes - `genomes_aligned.fasta`
* (optional, `--build-tree`) Generate a phylogenetic tree of the resolved genomes - `genomes_phylogenetic_tree.tree`

**Note:** Tree building is off by default -- it's slow at metagenome scale (thousands of
genomes), so most workflows should run phylogenetics as a separate downstream step on
`resolved_paths.fasta` instead. Pass `--build-tree` to enable it per-run. The tree file
`genomes_phylogenetic_tree.tree` is in newick format and can be visualised by tools such as
[iTOL (Interactive Tree of Life)](https://itol.embl.de/).


## Step-wise usage

You can execute each of the preprocessing, phables and postprocessing steps individually if you wish to do so as follows.

### Preprocessing only

You can use the following command to **only run the preprocessing steps**.

```bash
# Only preprocess data
phables run --input assembly_graph.gfa --reads fastq --threads 8 preprocess
```

### Genome resolution only

You can use the following command to **only run the genome resolution steps**. Please make sure to have the preprocessing results in the output folder.

```bash
# Only run phables core using short reads
phables run --input assembly_graph.gfa --reads fastq --threads 8 phables

# Only run phables core using long reads
phables run --input assembly_graph.gfa --reads fastq --threads 8 --longreads phables
```

### Postprocessing only

You can use the following command to **only run the postprocessing steps**.

```bash
# Only run phables core
phables run --input assembly_graph.gfa --reads fastq --threads 8 postprocess
```
