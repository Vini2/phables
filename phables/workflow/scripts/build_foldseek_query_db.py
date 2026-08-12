#!/usr/bin/env python3
"""
Build a foldseek query DB from an amino-acid FASTA + a ProstT5-predicted 3Di FASTA.

`foldseek createdb` on a plain protein FASTA fails ("No structures found in given
input") -- it expects real structure files. The correct approach, taken directly from
phold's own create_foldseek_db.generate_foldseek_db_from_aa_3di(), is `foldseek tsv2db`
on three numbered TSVs (AA, 3Di, header), matched by shared row position rather than by
building anything resembling real coordinates.

`tsv2db` on the AA sequences also auto-creates an empty `_ca` (C-alpha) companion file.
Its mere presence makes a later `foldseek search` try to run a TM-align step and fail
(``Structure alignment step died`` / ``getData: local id (4294967295) >= db size (0)``),
even though --alignment-type already defaults to 2 (3Di+AA, no coordinates needed). So
the empty `_ca` files are deleted here, immediately after the DB is built.
"""

import argparse
import subprocess
import sys
from pathlib import Path


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


def run(input_aa, input_3di, out_prefix, foldseek="foldseek"):
    aa = dict(read_fasta(input_aa))
    di = dict(read_fasta(input_3di))
    ids = [i for i in aa if i in di]
    missing = len(aa) - len(ids)
    if missing:
        print(
            f"WARNING: {missing} protein(s) have no 3Di prediction and will be "
            f"excluded from the query DB",
            file=sys.stderr,
        )
    if not ids:
        sys.exit("No proteins with both an amino-acid and a 3Di sequence -- nothing to build")

    out_prefix = Path(out_prefix)
    tmp_dir = out_prefix.parent
    aa_tsv = tmp_dir / f"{out_prefix.name}.aa.tsv"
    di_tsv = tmp_dir / f"{out_prefix.name}.3di.tsv"
    h_tsv = tmp_dir / f"{out_prefix.name}.header.tsv"

    with open(aa_tsv, "w") as f_aa, open(di_tsv, "w") as f_di, open(h_tsv, "w") as f_h:
        for i, seqid in enumerate(ids, 1):
            f_aa.write(f"{i}\t{aa[seqid]}\n")
            f_di.write(f"{i}\t{di[seqid]}\n")
            f_h.write(f"{i}\t{seqid}\n")

    def tsv2db(tsv_path, db_path, dbtype):
        subprocess.run(
            [foldseek, "tsv2db", str(tsv_path), str(db_path), "--output-dbtype", str(dbtype)],
            check=True,
        )

    tsv2db(aa_tsv, out_prefix, 0)
    tsv2db(di_tsv, f"{out_prefix}_ss", 0)
    tsv2db(h_tsv, f"{out_prefix}_h", 12)

    for f in aa_tsv, di_tsv, h_tsv:
        f.unlink()

    for suffix in ("", ".dbtype", ".index"):
        ca_file = Path(f"{out_prefix}_ca{suffix}")
        if ca_file.exists():
            ca_file.unlink()

    print(f"built query DB for {len(ids)} proteins -> {out_prefix}", file=sys.stderr)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aa", required=True, help="amino-acid protein FASTA")
    parser.add_argument("--threedi", required=True, help="ProstT5-predicted 3Di FASTA")
    parser.add_argument("--out-prefix", required=True, help="foldseek query DB prefix to create")
    parser.add_argument("--foldseek", default="foldseek")
    args = parser.parse_args(argv)
    run(args.aa, args.threedi, args.out_prefix, foldseek=args.foldseek)


if __name__ == "__main__":
    if "snakemake" in globals():
        run(
            snakemake.input.faa,  # noqa: F821
            snakemake.input.threedi,  # noqa: F821
            snakemake.params.out_prefix,  # noqa: F821
        )
    else:
        main()
