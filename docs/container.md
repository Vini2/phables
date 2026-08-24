# Running phables from a single container

`container/Dockerfile` builds **one monolithic image**: a ROCm base, conda
installed inside it, and every per-rule conda env this workflow can reach
already built into the image. Convert it to a single `.sif` and each
`phables run` inside it finds its environments already there.

There are no container-related CLI flags. Earlier versions had `--container`
and `--prostt5-container`, which pointed individual Snakemake rules at images
via `container:` directives; both are **removed**. The rules carry plain
`conda:` directives only, and the container satisfies them by having the envs
pre-built rather than by bypassing conda.

## Why

`--use-conda`'s per-rule env creation has a real failure mode at HPC scale.
Every rule's conda env is named by a hash of its env file, so every sample's
`phables run` needs the *same* env directory — but each sample is its own
separate Snakemake process (a SLURM array task), and Snakemake only serialises
env creation *within* one process's DAG. Two array tasks both needing the same
not-yet-built env at the same moment can race `mamba env create` into the same
directory and corrupt it (`Fatal Python error: init_fs_encoding ... no codec
search functions registered` is what that looks like in the wild).

Baking the envs into the image removes the failure mode entirely: at runtime
there is nothing left to create, so nothing left to race.

## Getting an image

CI (`.github/workflows/build_container.yaml`) builds and pushes to
`quay.io/<QUAY_USERNAME>/phables:<commit-sha>` on every push and pull request —
tag by the exact commit you want; there is no floating `latest`.

On Setonix, pull it once as a `.sif`:

```bash
module load singularity/4.1.0-slurm
singularity pull phables.sif docker://quay.io/gbouras13/phables:<commit-sha>
```

## Running it

Run phables *inside* the container. Nothing special is passed to phables
itself — `--use-conda` is its own default and the envs are already present:

```bash
singularity exec --rocm \
    -B /scratch/pawsey1018:/scratch/pawsey1018 \
    phables.sif \
    phables run --input assembly_graph.gfa --reads fastq \
        --output phables_out \
        --databases /scratch/.../all_databases/databases \
        --phagedetection prostt5-foldseek \
        --gpu-backend system \
        --prostt5-checkpoint /scratch/.../model.pt \
        --threads 8
```

Notes that matter on Setonix:

- **`--rocm`** exposes the host's GPU devices to the container. Without it,
  ProstT5 silently falls back to CPU.
- **Bind-mount your scratch** so databases, inputs and outputs are visible.
  Databases are deliberately *not* in the image (see below).
- **Don't pass `--conda-prefix`.** The envs were built at the default prefix
  (inside the installed package), and Snakemake resolves an env by hashing its
  file content *together with the prefix path* — changing the prefix changes
  the hash, and Snakemake would try to rebuild into a read-only filesystem.
- **Pass `--gpu-backend system` explicitly.** This is required, not optional.
  phables merges every CLI option over the config file
  (`merge_config=kwargs`), and `--gpu-backend`'s click default is `cpu` — so
  the image's own `config.yaml` value is *not* the effective value and cannot
  be relied on. Omitting the flag was tried and failed for real: the runtime
  config showed `gpu_backend: cpu`, Snakemake went to build a `prostt5-cpu`
  env that isn't in the image, and the run died with
  `OSError: [Errno 30] Read-only file system`. Any backend other than
  `system` fails the same way.

## What's in the image

- **Base**: `quay.io/pawsey/pytorch:2.7.1-rocm6.3.3`, Pawsey's own verified
  ROCm build, supplying both the ROCm userspace matching Setonix's MI250X
  (gfx90a) **and the torch the workflow actually uses**.
- **phables, Snakemake and pholdlib installed into that base python** — not
  into miniforge's. Snakemake runs a `script:` rule that declares no conda env
  using its own interpreter, so that interpreter has to be the one holding
  torch. The build asserts torch's version is unchanged across the pip install,
  since a silently-replaced torch is the exact failure this avoids.
- **Miniforge** at `/opt/miniforge3`, appended to `PATH` (never prepended, so
  it cannot shadow the base python). It exists only to provide the `conda`
  binary that builds the envs below.
- **Every per-rule conda env**, prebuilt by `container/prebuild_envs.sh`:
  coverm (minimap2/samtools/CoverM), genecall (FragGeneScan), pyrodigal-gv,
  smg (HMMER), mmseqs, foldseek, phylotree (MAFFT + cogent3/piqtree), and curl
  for `phables install`. **No `prostt5-*` env** — that's the point of
  `system`.

**Databases are not included** — PHROGs and the hallmark DB are multi-GB and
separately versioned. Mount them and point `--databases` at them, exactly as
outside a container. `phables install` also works inside the image if you'd
rather fetch them from there.

## Building it yourself

```bash
docker build -f container/Dockerfile -t phables:local .
```

`container/prebuild_envs.sh` does the env pre-building. It generates its own
throwaway inputs — a two-segment GFA and a pair of tiny gzipped FASTQs, plus
zero-byte placeholder database files — purely so the Snakemake DAG can
*resolve*, then runs `--conda-create-envs-only` once per flag combination (gene
caller, detection mode, GPU backend, tree) so every reachable env gets built.
Nothing is ever processed, since no job runs.

Two things it deliberately does **not** do, both of which broke a real build:

- It doesn't use `tests/data/`. That directory is in `.gitignore`, so it
  doesn't exist in a fresh clone or a CI checkout — depending on it failed with
  `Invalid value for '--reads': Path ... does not exist`.
- It points `phables install` at an *empty* databases directory, not the
  placeholder one. The placeholders satisfy install.smk's own download targets,
  so Snakemake reports "Nothing to be done", the DAG is empty and the `curl`
  env is silently never built.
`container/test_image.sh` then smoke-tests the result and fails the build if any
env is missing a binary the rules actually invoke, if torch/pholdlib aren't
importable in the ambient python, or if any pre-built env turns out to contain a
torch of its own — that last one being the regression that would mean a second
copy got installed after all.

It reports the base torch's build flavour (`torch.version.hip`) but does **not**
assert on it. The Pawsey base is a source build reporting e.g.
`2.7.1a0+gite2d141d`, with no `+rocm6.3` suffix; an earlier version of this
script tested for the substring `rocm` and failed a perfectly good ROCm build at
the very last step of a multi-GB image. What that torch is compiled against is
Pawsey's business.

**Disk**: this image is large — a ~14GB compressed ROCm base plus the conda
envs. The CI workflow runs `jlumbroso/free-disk-space` first because a stock
GitHub Actions runner has only ~14GB free and the base alone won't fit. Reusing
the base torch rather than installing a second one keeps several GB off the
total, but this is still close to the limit of what a hosted runner can build;
if CI starts failing on `no space left on device`, building on a machine with
real disk and pushing manually is the fallback.

## Status

A real `docker build` now gets all the way through the image and into
`test_image.sh`. Confirmed working on a real build:

- All **8** per-rule conda envs solve and install (coverm, genecall, smg,
  mmseqs, foldseek, phylotree, phables, curl) — including the pins that were
  only repodata-checked before (`mmseqs2=13.45111`, `cogent3<2026.7`).
- **No `prostt5-*` env is created**, so no second torch is downloaded — the
  torch reuse works as designed.
- Installing phables/Snakemake/pholdlib into the base python leaves its torch
  untouched (the build's own before/after assertion passed).
- Every per-rule binary resolves inside a pre-built env.

Still unverified: a Setonix Apptainer run — in particular whether ProstT5
actually sees the GPU through `singularity exec --rocm`, which no build-time
check can answer. Pull the `.sif`, put one real sample through it, and confirm
predict_3di lands on the GPU rather than silently falling back to CPU before
trusting this for production batches.
