#!/bin/bash
# Pre-builds EVERY per-rule conda env into the image, at Docker build time.
# Run from container/Dockerfile; not meant to be run on a host. Takes no
# arguments and reads nothing from the repo -- it generates its own throwaway
# inputs (see below), so it works in a bare CI checkout.
#
# Envs go into phables' own DEFAULT --conda-prefix (snake_base("workflow/conda"),
# i.e. inside the installed package) deliberately: a user inside the resulting
# container then runs plain `phables run ...` with no --conda-prefix flag and
# Snakemake resolves these exact envs. Snakemake names an env by hashing its
# file content together with the conda prefix path -- both identical at build
# time and run time here, since they're the same paths in the same image -- so
# the hashes match and nothing is rebuilt at runtime. That matters more than
# it sounds: a .sif is READ-ONLY when running, so an attempted rebuild is a
# hard failure, not just a slow path.
#
# --conda-create-envs-only builds a DAG's envs without running any of it. The
# DAG still has to RESOLVE, which requires the database files to EXIST -- but
# only to exist, since no job runs. Empty placeholder files are therefore
# enough. That was verified for real (dry-running every flag combination below
# against zero-byte placeholder DB files) before this script was written, and
# it's what keeps the multi-GB PHROGs/hallmark databases OUT of the image:
# mount the real ones at runtime via --databases, exactly as outside a
# container.
#
# One invocation per flag combination, because which envs a DAG needs depends
# on the flags -- gene caller, phage-detection mode and GPU backend each select
# different rules and env files. Together these cover every env under
# workflow/envs/ that any `phables run` (or `phables install`) can reach.

set -euxo pipefail

DB=/tmp/placeholder_db
WORK=/tmp/envbuild_inputs

mkdir -p "$DB/phrogs_mmseqs_db" "$DB/hallmark_db"
touch "$DB/marker.hmm" \
      "$DB/phrog_annot_v4.tsv" \
      "$DB/phrogs_mmseqs_db/phrogs_profile_db" \
      "$DB/hallmark_db/hallmark_db" \
      "$DB/hallmark_db/hallmark_categories.tsv"

# Synthetic minimal inputs, generated here rather than taken from
# tests/data/: that directory is in .gitignore, so it does NOT exist in a
# fresh clone or in a CI checkout -- depending on it made the image build
# fail with "Invalid value for '--reads': Path ... does not exist". Nothing
# below is ever actually processed (no job runs under
# --conda-create-envs-only); these files exist purely so click's exists=True
# checks pass and the DAG can resolve. Verified to produce the same DAG as
# the real test data.
mkdir -p "$WORK/reads"
printf 'H\tVN:Z:1.0\nS\tedge_1\tACGTACGTACGTACGTACGTACGTACGTACGT\tLN:i:32\nS\tedge_2\tTTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA\tLN:i:32\nL\tedge_1\t+\tedge_2\t+\t0M\n' \
    > "$WORK/assembly_graph.gfa"
printf '@r1\nACGT\n+\nIIII\n' | gzip > "$WORK/reads/sample1_R1.fastq.gz"
printf '@r1\nACGT\n+\nIIII\n' | gzip > "$WORK/reads/sample1_R2.fastq.gz"

GFA="$WORK/assembly_graph.gfa"
READS="$WORK/reads"
COMMON=(--input "$GFA" --reads "$READS" --databases "$DB" --threads 1)

# 1. default path -> coverm, genecall (FragGeneScan), smg (HMMER), mmseqs, phables
phables run "${COMMON[@]}" --output /tmp/envbuild1 --conda-create-envs-only

# 2. the other gene caller
phables run "${COMMON[@]}" --output /tmp/envbuild2 \
    --genecaller pyrodigal-gv --conda-create-envs-only

# 3. ProstT5 + foldseek detection -> the foldseek env (and, via
#    --gpu-backend system, NO prostt5-* torch env at all: predict_3di reuses the
#    base image's already-working ROCm torch, which is installed in the same
#    python running Snakemake here). This is the whole reason no multi-GB
#    second torch is downloaded during this build. Deliberately NOT `rocm`/
#    `cpu`/`cuda` -- each of those would solve and download their own torch.
phables run "${COMMON[@]}" --output /tmp/envbuild3 \
    --phagedetection prostt5-foldseek --gpu-backend system --conda-create-envs-only

# 4. Optional phylogenetic tree -> phylotree env (MAFFT + cogent3/piqtree).
#    Note conda resolves cogent3 fine, unlike pip, where every published
#    release is a prerelease and a plain version range matches nothing.
phables run "${COMMON[@]}" --output /tmp/envbuild4 \
    --build-tree --conda-create-envs-only

# 5. install.smk's own env (curl), so `phables install` works inside here too.
#    Deliberately pointed at an EMPTY databases dir, not "$DB": the placeholder
#    files in $DB satisfy install.smk's own download targets, so Snakemake says
#    "Nothing to be done", the DAG is empty, and NO env gets created -- the
#    curl env would silently be missing from the image. An empty dir puts the
#    four *_download rules in the DAG so their conda env actually gets built.
#    (--conda-create-envs-only still downloads nothing.)
mkdir -p /tmp/empty_db
phables install --output /tmp/envbuild5 --databases /tmp/empty_db --conda-create-envs-only

rm -rf /tmp/envbuild1 /tmp/envbuild2 /tmp/envbuild3 /tmp/envbuild4 /tmp/envbuild5 \
       "$DB" "$WORK" /tmp/empty_db
conda clean -a -y

echo "=== pre-built conda envs ==="
CONDA_PREFIX_DIR="$(python -c 'import phables, os; print(os.path.join(os.path.dirname(phables.__file__), "workflow", "conda"))')"
ls -1 "$CONDA_PREFIX_DIR"
