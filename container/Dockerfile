#
# phables -- ONE monolithic image for Setonix: ROCm base + conda + every
# per-rule conda env pre-built INSIDE the image. Build it once, convert it to a
# single .sif, and every `phables run` inside it finds its envs already there.
# Same strategy as hybracter's own container (miniforge installed into the
# image, envs materialised at build time), NOT the per-rule `container:`
# directive approach -- phables has no --container/--prostt5-container flags
# any more, and its Snakemake rules carry plain `conda:` directives only.
#
# Why this exists: a real, hit-in-production race condition. Every sample's
# `phables run` is its own separate process (one SLURM array task each), and
# Snakemake only serialises conda-env creation WITHIN a single process's DAG.
# Many array tasks all needing the same not-yet-built shared env at once
# corrupt each other's `mamba env create` ("Fatal Python error:
# init_fs_encoding ... no codec search functions registered"). With every env
# baked into the image there is nothing left to create at runtime, and so
# nothing left to race.
#
# TORCH: this image does NOT install a second PyTorch. The Pawsey ROCm base
# already ships a torch that is known to work on Setonix's MI250X (gfx90a), so
# predict_3di reuses it via `--gpu-backend system` (see below) instead of a
# prostt5-* conda env. That is the ONLY way to reuse it: conda envs are
# isolated, so a rule declaring `conda:` can never see a torch installed
# outside its env -- which is exactly why an earlier version of this image
# ended up downloading multiple GB of a second, unverified torch.
#
# The consequence for how the image is built: phables and Snakemake are
# installed into the BASE IMAGE'S python, not into miniforge's. Snakemake runs
# `script:` rules that declare no conda env using its own interpreter, so that
# interpreter has to be the one holding torch. Miniforge is installed only to
# provide the `conda` binary that BUILDS the other per-rule envs, and is
# appended to PATH rather than prepended so it can never shadow the base
# python.
#
# See ../docs/container.md for how to build the .sif and run it on Setonix.
#
FROM quay.io/pawsey/pytorch:2.7.1-rocm6.3.3

ARG DEBIAN_FRONTEND="noninteractive"

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget \
        bzip2 \
        ca-certificates \
        git \
        procps \
    && rm -rf /var/lib/apt/lists/*

# Pin down the base image's python BEFORE miniforge exists, so later steps can
# target it unambiguously no matter where the base keeps it. Everything phables
# needs goes into this interpreter.
RUN command -v python3 > /opt/base_python_path && \
    ln -sf "$(cat /opt/base_python_path)" /usr/local/bin/python && \
    echo "base python: $(cat /opt/base_python_path)" && \
    "$(cat /opt/base_python_path)" -c "import torch; print('base torch:', torch.__version__)"

#
# Miniforge -- installed ONLY for its `conda`/`mamba` binaries, which Snakemake
# shells out to when building the per-rule envs. Deliberately APPENDED to PATH:
# prepending it would make miniforge's python the default, and phables would
# then run under an interpreter with no torch in it, silently defeating the
# whole point of building on a ROCm base.
#
ENV CONDA_DIR=/opt/miniforge3
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p "${CONDA_DIR}" && \
    rm /tmp/miniforge.sh
ENV PATH="${PATH}:${CONDA_DIR}/bin"

#
# phables + Snakemake + pholdlib into the BASE python. pholdlib declares no
# torch version constraint, so pip pulls its own dependencies (transformers,
# sentencepiece, h5py, ...) without touching the already-installed torch -- the
# assertion below fails the build if that ever stops being true, since a
# silently-replaced torch is precisely the failure this image exists to avoid.
#
COPY . /opt/phables_src
RUN BASE_PY="$(cat /opt/base_python_path)" && \
    TORCH_BEFORE="$("$BASE_PY" -c 'import torch; print(torch.__version__)')" && \
    "$BASE_PY" -m pip install --no-cache-dir "snakemake>=8" pholdlib /opt/phables_src && \
    TORCH_AFTER="$("$BASE_PY" -c 'import torch; print(torch.__version__)')" && \
    echo "torch before=${TORCH_BEFORE} after=${TORCH_AFTER}" && \
    test "$TORCH_BEFORE" = "$TORCH_AFTER"

#
# NOTE: there is deliberately NO step here patching config.yaml's gpu_backend
# to `system`. An earlier version did exactly that, and it does not work:
# phables passes every CLI option to snaketool as merge_config=kwargs, so any
# option with a non-None click default is merged OVER the config file, and
# --gpu-backend's click default is "cpu". The patched config was silently
# ignored, Snakemake went to build a prostt5-cpu env that isn't in this image,
# and the run died on the read-only .sif filesystem.
#
# Callers must therefore pass `--gpu-backend system` explicitly. See
# ../docs/container.md; mass_processing/run_phables_container.sh does this.
#

# Pre-build every per-rule conda env. See the script for the full rationale
# (default --conda-prefix, placeholder databases, one pass per flag combo).
RUN bash /opt/phables_src/container/prebuild_envs.sh

# Fail the build loudly rather than pushing a broken image to quay.io.
RUN bash /opt/phables_src/container/test_image.sh
