#!/usr/bin/env python3
"""
ProstT5 -> 3Di prediction for unitig ORF proteins.

Reuses pholdlib (github.com/gbouras13/phold-lib) directly rather than reimplementing
ProstT5 batching/fp16 handling -- it is the shared inference engine phold itself wraps.

Deliberately does NOT import phold: phold's get_embeddings() expects a nested
{contig_id: {seq_id: BioPython_feature}} cds_dict plus a GenBank-oriented database
layer for weight downloads (phold.databases.db). That machinery is unneeded for
"predict 3Di for a flat protein FASTA" -- this script calls pholdlib's own
get_T5_model / load_predictor / run_prostt5_inference directly.

Input protein headers are expected in the {unitig}_{start}_{end}_{strand} convention
produced by gene_caller.py, so the output 3Di FASTA carries the same headers straight
through to `foldseek search` and downstream unitig-name recovery (hallmark_utils.py)
stays consistent with the marker-gene path.
"""

import argparse
import sys
from pathlib import Path

# Work around transformers' loss registry unconditionally importing torchaudio.
# Root cause: merely importing pholdlib.prostt5 executes pholdlib/prostt5/model.py,
# which does `from transformers import T5EncoderModel, ...` -- and recent
# transformers versions route model imports through a shared loss registry
# (transformers/loss/loss_utils.py) that, once it decides torchaudio is
# "available", imports transformers/loss/loss_rnnt.py to support
# ParakeetForRNNTLoss -- an audio ASR loss completely unrelated to T5/protein
# embeddings -- which does a plain `import torchaudio`. We never touch anything
# audio-related, so if that import can't actually succeed, install a harmless
# stub in sys.modules so transformers' import chain completes anyway. Two
# distinct real failure modes hit in practice, both caught here:
#   - OSError: torchaudio is installed but its compiled extension can't load
#     (confirmed on Setonix's ROCm container: `libomp.so: cannot open shared
#     object file` -- the image doesn't have that runtime linkable).
#   - ModuleNotFoundError (a subclass of ImportError): torchaudio isn't
#     installed at all, e.g. a plain conda env that only installs torch +
#     pholdlib, not the full audio/vision extras a generic "pytorch" container
#     image tends to bundle. Confirmed this needed covering too, not just
#     OSError -- our own diagnostic `import torchaudio` below would otherwise
#     crash unnecessarily in exactly this env, independent of whether
#     transformers' own internal logic would have needed torchaudio at all.
# No-op everywhere torchaudio genuinely works normally.
#
# The stub needs a real (if empty) __spec__, not just a bare ModuleType: before
# transformers ever reaches the `import torchaudio` line above, it first runs a
# *lightweight* availability check (is_torchaudio_available() ->
# importlib.util.find_spec("torchaudio")) that explicitly raises
# `ValueError: torchaudio.__spec__ is None` if sys.modules already has an entry
# for the name with no spec -- confirmed by hitting this exact error with a first
# version of this workaround that used a bare ModuleType. __version__ is set too,
# defensively: if the real package's installed-distribution metadata ever isn't
# resolvable, transformers falls back to reading torchaudio.__version__ directly,
# which our stub would otherwise lack.
try:
    import torchaudio  # noqa: F401
except (OSError, ImportError):
    import importlib.machinery
    import types

    _torchaudio_stub = types.ModuleType("torchaudio")
    _torchaudio_stub.__spec__ = importlib.machinery.ModuleSpec(
        "torchaudio", loader=None
    )
    _torchaudio_stub.__version__ = "0.0.0"
    sys.modules["torchaudio"] = _torchaudio_stub

from pholdlib.prostt5.device import parse_gpus
from pholdlib.prostt5.model import get_T5_model, load_predictor
from pholdlib.prostt5.inference import run_prostt5_inference
from pholdlib.prostt5.output import SS_MAPPING


def read_fasta(path):
    header, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def hf_cache_present(model_dir, model_name):
    """True if the HF snapshot already looks downloaded -- avoids forcing
    local_files_only=True (pholdlib's default with no check_fn) on a genuine first
    run, which raises rather than downloading. See get_T5_model in pholdlib.prostt5.model.
    """
    cache_path = Path(model_dir) / f"models--{model_name.replace('/', '--')}"
    snapshots = cache_path / "snapshots"
    if not snapshots.exists():
        return False
    return any(snapshots.iterdir())


def run(
    input_fasta,
    output_3di,
    checkpoint,
    model_name,
    model_dir,
    threads=1,
    half_precision=False,
    cpu=False,
    max_residues=100000,
    max_seq_len=30000,
    max_batch=10000,
    mask_threshold=0.0,
):
    if not model_dir:
        # Snakemake's script: entrypoint passes config["prostt5_model_dir"] straight
        # through, which is None unless --prostt5-model-dir was explicitly set --
        # unlike the CLI entrypoint's own --model-dir, which argparse defaults to
        # this same path. Both entrypoints need the fallback, so it lives here.
        model_dir = str(Path.home() / ".cache" / "prostt5")

    records = list(read_fasta(input_fasta))
    if not records:
        sys.exit(f"No sequences read from {input_fasta}")
    seq_dict = [(rid, seq, len(seq)) for rid, seq in records]
    seq_dict.sort(key=lambda x: x[2], reverse=True)

    need_download = not hf_cache_present(model_dir, model_name)
    devices = parse_gpus(cpu, None)
    device_str = devices[0]
    print(
        f"device: {device_str}  |  {len(seq_dict)} proteins  |  "
        f"download needed: {need_download}",
        file=sys.stderr,
    )

    model, vocab, device = get_T5_model(
        model_dir=model_dir,
        model_name=model_name,
        cpu=cpu,
        threads=threads,
        check_fn=(lambda *_: need_download),
        device=None if device_str in ("cuda:0",) else device_str,
    )

    if half_precision:
        if device.type == "cpu":
            print("CPU device -- ignoring half_precision", file=sys.stderr)
        else:
            model = model.half()

    predictor = load_predictor(checkpoint, device)
    if half_precision and device.type != "cpu":
        predictor = predictor.half()

    preds, _, _, fail_ids = run_prostt5_inference(
        seq_dict,
        model,
        vocab,
        predictor,
        device,
        max_residues=max_residues,
        max_seq_len=max_seq_len,
        max_batch=max_batch,
        output_probs=True,
        desc="Predicting 3Di",
    )

    if fail_ids:
        print(f"WARNING: {len(fail_ids)} sequences failed: {fail_ids[:5]}", file=sys.stderr)

    mask_prop = mask_threshold / 100
    n_written = 0
    with open(output_3di, "w") as out:
        for rid, _ in records:
            if rid not in preds:
                continue
            pred, mean_prob, all_prob = preds[rid]
            if len(pred) == 0:
                continue
            pred = pred.copy()
            pred[all_prob[0] < mask_prop] = 20  # 'X'
            threeDi = "".join(SS_MAPPING[int(y)] for y in pred)
            out.write(f">{rid}\n{threeDi}\n")
            n_written += 1

    print(f"wrote {n_written} 3Di sequences -> {output_3di}", file=sys.stderr)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, help="protein FASTA")
    parser.add_argument("-o", "--output", required=True, help="output 3Di FASTA")
    parser.add_argument("--checkpoint", required=True, help="CNN prediction-head .pt/.pth")
    parser.add_argument("--model-name", default="Rostlab/ProstT5_fp16")
    parser.add_argument("--model-dir", default=str(Path.home() / ".cache" / "prostt5"))
    parser.add_argument("--half-precision", action="store_true")
    parser.add_argument("--cpu", action="store_true")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--max-residues", type=int, default=100000)
    parser.add_argument("--max-seq-len", type=int, default=30000)
    parser.add_argument("--max-batch", type=int, default=10000)
    parser.add_argument("--mask-threshold", type=float, default=0.0)
    args = parser.parse_args(argv)
    run(
        args.input,
        args.output,
        args.checkpoint,
        args.model_name,
        args.model_dir,
        threads=args.threads,
        half_precision=args.half_precision,
        cpu=args.cpu,
        max_residues=args.max_residues,
        max_seq_len=args.max_seq_len,
        max_batch=args.max_batch,
        mask_threshold=args.mask_threshold,
    )


if __name__ == "__main__":
    # Snakemake's script: directive injects a `snakemake` object into globals rather
    # than passing argv, so support both entry points.
    if "snakemake" in globals():
        run(
            snakemake.input.faa,  # noqa: F821
            snakemake.output.threedi,  # noqa: F821
            snakemake.params.checkpoint,  # noqa: F821
            snakemake.params.model_name,  # noqa: F821
            snakemake.params.model_dir,  # noqa: F821
            threads=snakemake.threads,  # noqa: F821
            half_precision=snakemake.params.half_precision,  # noqa: F821
            cpu=snakemake.params.cpu,  # noqa: F821
            max_residues=snakemake.params.max_residues,  # noqa: F821
            max_seq_len=snakemake.params.max_seq_len,  # noqa: F821
            max_batch=snakemake.params.max_batch,  # noqa: F821
        )
    else:
        main()
