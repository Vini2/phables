#!/usr/bin/env python3
"""
Call genes on unitigs with pyrodigal-gv, as a drop-in replacement for FragGeneScan.

Protein headers deliberately reproduce FragGeneScan's `{seqid}_{start}_{end}_{strand}`
convention. phables_utils.gene_utils.get_smg_unitigs recovers the unitig name by
stripping the last three underscore-separated fields, so any other naming scheme
silently resolves every unitig to the empty string and disables marker-gene
filtering entirely. Keeping the convention lets both callers share one parser and
makes --gene-caller a genuine one-variable ablation.
"""

import argparse
import multiprocessing.pool
import sys

import pyrodigal_gv


def read_fasta(path):
    """Yield (id, sequence) pairs. The id is the header up to the first whitespace."""
    header = None
    chunks = []
    with open(path, "r") as handle:
        for line in handle:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def gene_id(unitig, gene):
    """Build a FragGeneScan-compatible protein id."""
    strand = "+" if gene.strand == 1 else "-"
    return f"{unitig}_{gene.begin}_{gene.end}_{strand}"


def call_genes(records, threads, viral_only=False):
    """Predict genes for (id, sequence) records. Returns (unitig, Genes) pairs."""
    # closed=False keeps genes that run off the end of a unitig. Unitigs are graph
    # fragments and truncate genes constantly; FragGeneScan was run with -complete=0
    # for the same reason.
    finder = pyrodigal_gv.ViralGeneFinder(
        meta=True,
        closed=False,
        viral_only=viral_only,
    )

    if threads > 1:
        with multiprocessing.pool.ThreadPool(threads) as pool:
            genes = pool.map(finder.find_genes, [seq for _, seq in records])
    else:
        genes = [finder.find_genes(seq) for _, seq in records]

    return list(zip([name for name, _ in records], genes))


def write_proteins(results, path):
    n = 0
    with open(path, "w") as handle:
        for unitig, genes in results:
            for gene in genes:
                # include_stop=False matches FragGeneScan, which does not emit the
                # trailing stop codon.
                handle.write(f">{gene_id(unitig, gene)}\n")
                handle.write(f"{gene.translate(include_stop=False)}\n")
                n += 1
    return n


def run(input_fasta, output_faa, threads, viral_only=False):
    records = list(read_fasta(input_fasta))
    if not records:
        sys.exit(f"No sequences read from {input_fasta}")

    results = call_genes(records, threads, viral_only=viral_only)
    n = write_proteins(results, output_faa)

    print(f"Called {n} genes across {len(records)} unitigs -> {output_faa}", file=sys.stderr)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, help="unitig FASTA (edges.fasta)")
    parser.add_argument("-o", "--output", required=True, help="output protein FASTA")
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument(
        "--viral-only",
        action="store_true",
        help=(
            "restrict to viral models. Off by default: the same proteins feed the "
            "bacterial single-copy marker gene search, which needs chromosomal genes called."
        ),
    )
    args = parser.parse_args(argv)
    run(args.input, args.output, args.threads, viral_only=args.viral_only)


if __name__ == "__main__":
    # Snakemake's script: directive injects a `snakemake` object into globals rather
    # than passing argv, so support both entry points.
    if "snakemake" in globals():
        run(
            snakemake.input.genome,  # noqa: F821
            snakemake.output.faa,  # noqa: F821
            snakemake.threads,  # noqa: F821
        )
    else:
        main()
