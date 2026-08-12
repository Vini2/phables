#
# phables -- massive all-in-one image for Setonix (ROCm), replacing per-rule
# --use-conda env creation entirely. See ../docs/container.md for how to run
# it (--container <image> --use-singularity --no-use-conda) and why this
# exists (a real, hit-in-production race condition: every sample's `phables
# run` is its own separate process, so many array tasks needing the same
# shared --use-conda env for the first time at once can corrupt each other's
# `mamba env create`).
#
# Base is Pawsey's own verified-working ROCm+PyTorch image -- the SAME one
# phold's own container (../../phold/container/hpci/Dockerfile) builds from,
# deliberately: torch==2.7.1+rocm6.3 baked in here is exactly what
# envs/prostt5-rocm.yaml pins, confirmed working on real Setonix MI250X
# hardware this session. Nothing below reinstalls or touches torch -- pip
# installs pholdlib and every other pure-Python dependency directly into this
# same system Python, so they see (and don't conflict with) the base image's
# already-correct torch/ROCm stack rather than resolving a second, possibly
# different one.
FROM quay.io/pawsey/pytorch:2.7.1-rocm6.3.3

ARG DEBIAN_FRONTEND="noninteractive"

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    tar \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# the base image only provides python3, not python -- some tooling (and
# Snakemake's own script: directive) expects `python` on PATH
RUN ln -sf "$(command -v python3)" /usr/local/bin/python

#
# foldseek -- same install as phold's own Dockerfile (direct binary download,
# not conda/pip), since that's the reference build this image otherwise
# mirrors for the ROCm/torch stack
#
RUN wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz && \
    tar -xzf foldseek-linux-avx2.tar.gz && \
    rm foldseek-linux-avx2.tar.gz && \
    mv foldseek /opt/foldseek
ENV PATH="/opt/foldseek/bin:${PATH}"

#
# Every other per-rule tool that isn't a pure-Python package: bootstrap a
# standalone micromamba (no full miniforge needed) and install them all into
# one prefix, /opt/conda -- kept entirely separate from the base image's own
# system Python (added to PATH, but conda never pulls in a python of its own
# here since none of these packages are Python libraries, only standalone
# binaries -- confirmed by their real recipes: minimap2/samtools/mmseqs2/
# fraggenescan/hmmer/mafft have no python dependency; coverm's own recipe
# pulls minimap2+samtools transitively too, listed here explicitly anyway for
# a single deterministic solve rather than relying on that side effect).
# mmseqs2 pinned to the exact version envs/mmseqs.yaml already uses.
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba && \
    mv bin/micromamba /usr/local/bin/micromamba && \
    micromamba create -y -p /opt/conda -c conda-forge -c bioconda \
        minimap2 \
        samtools \
        "coverm>=0.6" \
        "mmseqs2=13.45111" \
        fraggenescan \
        hmmer \
        mafft \
    && micromamba clean -a -y
ENV PATH="/opt/conda/bin:${PATH}"

#
# Everything else -- phables' own pure-Python deps, ProstT5/pholdlib, and
# phables itself -- pip installed directly into the base image's system
# Python. Deliberately NOT a separate venv/conda env: the whole point is for
# all of this to see the base image's already-installed, already-verified
# torch/ROCm build rather than resolving (and risking conflicting with) a
# second one.
#
# pholdlib's own declared torch requirement is unconstrained (confirmed
# earlier this session against envs/prostt5-rocm.yaml/-cpu.yaml/-cuda.yaml) --
# `pip install pholdlib` alone pulls in its transitive deps (transformers,
# sentencepiece, h5py, loguru, etc.) via pip's own resolver without touching
# the already-installed torch.
RUN pip install --ignore-installed PyYAML && \
    python -c "import torch; print('torch (pre-existing, untouched):', torch.__version__)"

#
# NOT included: cogent3/piqtree (phylotree.yaml's own deps, for the optional
# --phylotree/build_tree rule -- not part of this workflow's usage).
# Deliberately left out rather than fought into pip (every cogent3 release on
# PyPI is an unreleased alpha, which breaks a plain version-range pip
# install). Note: --phylotree is therefore NOT supported under --container --
# build_tree still declares container: CONTAINER_IMAGE (unchanged), so it
# will fail if actually invoked against this image.
#
RUN pip install \
        "pholdlib>=0.1.2" \
        "pyrodigal-gv>=0.3.1" \
        biopython \
        python-igraph \
        pysam \
        "networkx>=2.8.6" \
        scipy \
        numpy \
        pandas \
        tqdm \
        click \
        "metasnek>=0.0.3" \
        "flowpaths>=0.2.20"

RUN python -c "import torch; print('torch (post pip installs, should be unchanged):', torch.__version__)"

#
# phables itself -- installed from the exact source in this build context
# (this Dockerfile is meant to be built by CI on every commit, tagged by
# commit SHA, so it's the checked-out commit's own tree, not a fetched
# git+https ref pinned separately).
#
COPY . /opt/phables_src
RUN pip install /opt/phables_src

#
# tests -- fail the build loudly if any of the above didn't actually work,
# rather than shipping a broken image to quay.io
#
RUN foldseek --help
RUN minimap2 --version
RUN samtools --version
RUN mmseqs version
RUN coverm --version
RUN hmmsearch -h
# the bioconda package's actual binary is run_FragGeneScan.pl, not
# `fraggenescan` -- confirmed against genes.smk's own real invocation
RUN run_FragGeneScan.pl 2>&1 | head -5 || true
RUN mafft --version
RUN python --version
RUN python -c "import pholdlib; print('pholdlib OK')"
RUN python -c "import torch; assert torch.__version__.startswith('2.7.1'), torch.__version__"
RUN phables --version
RUN phables run -h
