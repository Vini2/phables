#!/bin/bash
# Build-time smoke tests for the monolithic image -- fail the build loudly
# rather than pushing something broken to quay.io. Run from container/Dockerfile.
#
# Every per-rule tool is checked INSIDE whichever pre-built conda env provides
# it, not on the base PATH: none of them are on the base PATH, and that's the
# point of the per-rule env layout. Envs are located by searching the prefix
# rather than by hardcoded directory names, because Snakemake names each env
# by a content hash that this script has no business predicting.

set -euo pipefail

PREFIX="$(python -c 'import phables, os; print(os.path.join(os.path.dirname(phables.__file__), "workflow", "conda"))')"
echo "conda prefix: $PREFIX"
test -d "$PREFIX"

echo "=== phables CLI ==="
phables --version
phables run -h > /dev/null
phables install -h > /dev/null
echo "CLI OK"

echo "=== the --container / --prostt5-container flags must be GONE ==="
if phables run -h 2>&1 | grep -qE '\-\-(prostt5-)?container'; then
    echo "ERROR: a container flag is still present in the CLI" >&2
    exit 1
fi
echo "confirmed absent"

echo "=== pre-built env count ==="
n=$(find "$PREFIX" -maxdepth 1 -mindepth 1 -type d | wc -l)
echo "found $n env directories"
find "$PREFIX" -maxdepth 1 -mindepth 1 -type d -exec basename {} \;
# 6 distinct env files are reachable at minimum: coverm, genecall, smg, mmseqs,
# phables, curl -- before counting foldseek and phylotree. No prostt5-* env is
# expected (or wanted): gpu_backend=system reuses the base image's torch.
test "$n" -ge 6

# Finds an executable in any pre-built env; fails if no env provides it.
check_bin() {
    local want="$1" e
    for e in "$PREFIX"/*/; do
        if [ -x "${e}bin/${want}" ]; then
            echo "OK: $want -> $e"
            return 0
        fi
    done
    echo "MISSING from every pre-built env: $want" >&2
    return 1
}

echo "=== per-rule binaries ==="
check_bin minimap2
check_bin samtools
check_bin coverm
check_bin mmseqs
check_bin foldseek
check_bin hmmsearch
# the bioconda package's real binary name, confirmed against genes.smk's own
# invocation -- NOT `fraggenescan`
check_bin run_FragGeneScan.pl
check_bin mafft
check_bin curl

echo "=== torch must come from the BASE image, reused -- not reinstalled ==="
# predict_3di runs with NO conda env (gpu_backend=system), i.e. in the same
# python that runs Snakemake -- which must therefore be the base image's python,
# the one already holding a working torch. Two things are checked: that torch +
# pholdlib are importable here (fatal if not -- predict_3di simply cannot run),
# and that no conda env carries a torch of its own (fatal -- that would mean a
# second copy got installed after all). What the torch's build flavour is, is
# only reported.
#
# Reported, NOT asserted. This deliberately cannot fail the build.
#
# The base image's torch is Pawsey's own source build for Setonix -- it reports
# e.g. "2.7.1a0+gite2d141d", with no "+rocm6.3" suffix, because it isn't a
# stock wheel. An earlier version of this script tested for the substring
# "rocm" in torch.__version__ and failed a perfectly good ROCm build at the
# very last step of a multi-GB image. Whether that torch is ROCm-enabled is
# Pawsey's business, not something worth re-litigating here at build time, so
# this prints the facts (torch.version.hip is the real signal: set to the HIP
# version on ROCm builds, None otherwise) and moves on.
#
# The check that DOES matter -- that nothing installed a second torch -- is
# below and is fatal.
python - <<'PY'
import torch

hip = getattr(torch.version, "hip", None)
print("torch:", torch.__version__)
print("torch.version.hip:", hip)
print("torch.version.cuda:", getattr(torch.version, "cuda", None))
if hip is None:
    print("WARNING: torch.version.hip is None -- this torch does not look "
          "ROCm-enabled. Fine if intentional (e.g. a CPU-only base image); "
          "worth a look if you expected GPU ProstT5 on Setonix.")
PY
python -c "import pholdlib; print('pholdlib OK')"
python -c "import phables, snakemake; print('phables + snakemake share this interpreter')"

# A torch inside any pre-built env means a second copy got installed after all,
# which is the exact regression this layout exists to prevent.
for e in "$PREFIX"/*/; do
    if [ -x "${e}bin/python" ] && "${e}bin/python" -c "import torch" 2>/dev/null; then
        echo "ERROR: a pre-built conda env contains its own torch: $e" >&2
        echo "       gpu_backend=system should mean no prostt5-* env is built." >&2
        exit 1
    fi
done
echo "OK: no pre-built env ships a duplicate torch"

echo "=== --gpu-backend system must be available ==="
# This, NOT a grep of config.yaml's gpu_backend value. An earlier version of
# this script asserted the config file said "system" and passed happily, while
# actual runs still used cpu: phables merges every CLI option over the config
# (merge_config=kwargs), and --gpu-backend's click default is "cpu". The config
# file value is simply not the effective value, so checking it proves nothing.
# What the image can meaningfully guarantee is that the choice EXISTS -- the
# caller is responsible for passing it (see docs/container.md).
phables run -h 2>&1 | grep -q -- '--gpu-backend .*system' \
    || { echo "ERROR: this phables has no --gpu-backend system choice" >&2; exit 1; }
echo "OK: --gpu-backend system is accepted"

echo "=== all image tests passed ==="
