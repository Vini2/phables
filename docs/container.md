# Running phables from a single container

`--use-conda`'s per-rule environment creation is convenient but has a real
failure mode at HPC scale: every rule's conda env is named by a hash of its
env file, so every sample's `phables run` needs the *same* env directory —
but each sample is its own separate Snakemake process (a SLURM array task,
say), and Snakemake only serializes env creation *within* one process's DAG.
Two array tasks both needing the same not-yet-built env at the same moment
can race `mamba env create` into the same directory, corrupting it
(`Fatal Python error: init_fs_encoding ... no codec search functions
registered` is what that looks like when it happens).

`container/Dockerfile` builds a single image with every per-rule tool this
workflow needs already installed — `--container` then points the whole
workflow at it, replacing conda entirely so there's nothing left to race.

## Getting an image

CI (`.github/workflows/build_container.yaml`) builds and pushes an image to
`quay.io/<QUAY_USERNAME>/phables:<commit-sha>` on every push and pull request
— tag by the exact commit you want, there's no floating `latest`.

## Running with it

```bash
phables run --input assembly_graph.gfa --reads fastq \
    --container quay.io/gbouras13/phables:<commit-sha> \
    --no-use-conda \
    --use-singularity
```

`--use-singularity` is Snakemake's own flag for running rules inside
Apptainer/Singularity (the standard container runtime on HPC, incl. Setonix)
— pass it through phables' existing `snake_args` passthrough, same as
`--dry-run`/`--keep-going`/etc. `--use-conda` alone does **not** activate
container execution; explicitly turn it off (`--no-use-conda`) when using
`--container`, since combining the two makes Snakemake build a *separate*
conda env *inside* the container instead of using what's already
installed there — exactly the redundant-torch-reinstall problem this image
exists to avoid.

## What's actually in the image

Base: `quay.io/pawsey/pytorch:2.7.1-rocm6.3.3` — the same verified-working
ROCm+PyTorch build `envs/prostt5-rocm.yaml` pins, and the same base phold's
own container (`../../phold/container/hpci/Dockerfile`) builds from. Nothing
in phables' own Dockerfile reinstalls or touches torch — `pip install
pholdlib` and everything else pure-Python go straight into that same system
Python, so they use the base image's already-correct torch/ROCm stack rather
than resolving a second, possibly conflicting one. This is the direct answer
to "don't rebuild the ProstT5 conda env for this container" — there's no
separate ProstT5 env in the container at all; predict_3di runs against the
same Python everything else does.

Everything else — foldseek (direct binary, matching phold's own install),
minimap2/samtools/coverm/mmseqs2/FragGeneScan/HMMER/MAFFT (via a standalone
micromamba install into `/opt/conda`, kept deliberately separate from the
system Python so it can't touch it), and phables' own Python dependencies
(pyrodigal-gv, biopython, python-igraph, pysam, flowpaths, cogent3/piqtree,
...) — covers every `envs/*.yaml` this workflow has, so the whole DAG can run
from this one image with `--no-use-conda`.

## `--container` vs `--prostt5-container`

`--prostt5-container` already existed for pointing *just* `predict_3di` at a
container (e.g. phold's own image, which also has pholdlib+torch). It still
works exactly as before, and if set, it wins for that one rule. `--container`
is new and broader: it's the default container for **every** rule, including
predict_3di if `--prostt5-container` isn't also given. In practice, setting
`--container` alone is enough — there's no reason to use `--prostt5-container`
separately unless you specifically want predict_3di on a *different* image
than the rest of the workflow.

## Building it yourself

```bash
docker build -f container/Dockerfile -t phables:local .
```

The base image is large (~14GB compressed) — building on a standard GitHub
Actions runner needs the disk-cleanup step already in
`build_container.yaml` (`jlumbroso/free-disk-space`), or the very first
`FROM` line fails with `no space left on device`. Building locally needs
proportionate free disk space too.

**Not yet built or run for real** — the Dockerfile and the workflow-wide
`--container` plumbing were written and the Snakemake DAG wiring was verified
via real `--dry-run`s (with and without `--container` set, confirming both
resolve identically at the structural level), but no actual `docker build`
or real Apptainer/Singularity execution against Setonix hardware has
happened yet. Build it, push a first tag, and run one real sample through it
before trusting this for production batches.
