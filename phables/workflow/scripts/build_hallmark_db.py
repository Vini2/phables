#!/usr/bin/env python3
"""
Build the phage-hallmark Foldseek subDB that `--phagedetection prostt5-foldseek`
searches against (--hallmark-db / --hallmark-categories).

See docs/hallmark_db.md for the full explanation of why these categories, why
the representative cap, and the real numbers from the reference build this
script produces. Short version: subsets the full phold structure DB
(all_phold_structures, ~1.36M entries, github.com/gbouras13/phold) down to the
structural/virion PHROG categories -- head and packaging, connector, tail,
lysis -- which is what phables_utils/component_utils.py's get_components
actually screens for. integration and excision is a strong lysogeny signal but
integrases/recombinases are abundant on bacterial chromosomes and MGEs too, so
it's kept as a separate id list/subDB, not folded into the hallmark set.

phold_annots.tsv uses CRLF line endings; every field comparison below strips a
trailing \\r -- found the hard way (a first pass silently matched zero rows).

Usage:
    build_hallmark_db.py --phold-db-prefix all_phold_structures \\
        --annots phold_annots.tsv --out-dir hallmark_db/ --max-per-phrog 30
"""

import argparse
import random
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

__author__ = "George Bouras"
__copyright__ = "Copyright 2026, Phables Project"
__license__ = "MIT"
__type__ = "Support Script"
__maintainer__ = "George Bouras"
__email__ = "george.bouras@adelaide.edu.au"

HALLMARK_CATEGORIES = {"head and packaging", "connector", "tail", "lysis"}
SEPARATE_CHANNELS = {"integration and excision"}


def read_annots(path):
    """Yield (phrog_id, product, function) with CRLF stripped. Non-numeric phrog ids
    (acr, card, vfdb, defensefinder, netflax, ...) belong to other phold DB channels,
    not the PHROG catalogue, and are skipped."""
    with open(path, "r", newline="") as handle:
        header = handle.readline()
        assert header.startswith("phrog\t"), f"unexpected header: {header!r}"
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                continue
            phrog_id, product, function = fields[0], fields[1], fields[2].strip('"')
            if not phrog_id.isdigit():
                continue
            yield phrog_id, product, function


def categorise(annots_path):
    hallmark = {}
    separate = {}
    for phrog_id, product, function in read_annots(annots_path):
        if function in HALLMARK_CATEGORIES:
            hallmark[phrog_id] = function
        elif function in SEPARATE_CHANNELS:
            separate[phrog_id] = function
    return hallmark, separate


def load_lookup_by_phrog(lookup_path, wanted_ids):
    """Group lookup row indices by phrog id, for exactly the ids we need.

    Format: <row_index>\\t<prefix_phrogid:name>\\t<raw_index>. Only phrog_-prefixed
    entries are PHROG catalogue members (other prefixes are separate phold DB channels:
    acr, card, vfdb, defensefinder, dgr_phrog, efam_phrog, envhog_phrog, netflax).
    """
    by_phrog = defaultdict(list)
    with open(lookup_path, "r") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[1].startswith("phrog_"):
                continue
            key = fields[1].split(":", 1)[0]  # "phrog_<id>"
            phrog_id = key[len("phrog_"):]
            if phrog_id in wanted_ids:
                by_phrog[phrog_id].append(int(fields[0]))
    return by_phrog


def cap_representatives(by_phrog, max_per_phrog, seed):
    """Cap each PHROG family at max_per_phrog members.

    A deterministic stride sample rather than the first N: lookup order tends to cluster
    by source dataset, so a stride spreads the cap across whatever diversity exists in the
    family instead of taking N near-duplicates from one source.
    """
    rng = random.Random(seed)
    capped = {}
    for phrog_id, indices in by_phrog.items():
        if len(indices) <= max_per_phrog:
            capped[phrog_id] = list(indices)
            continue
        ordered = sorted(indices)
        stride = len(ordered) / max_per_phrog
        picked = [ordered[int(i * stride)] for i in range(max_per_phrog)]
        capped[phrog_id] = picked
    return capped


def run(cmd):
    print(f"$ {' '.join(cmd)}", file=sys.stderr)
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--phold-db-prefix",
        required=True,
        help="prefix of the full phold structure DB, e.g. all_phold_structures "
        "(expects <prefix>, <prefix>_ss, <prefix>_h, and <prefix>.lookup alongside it)",
    )
    parser.add_argument(
        "--annots", required=True, help="phold_annots.tsv (bundled with phold's own DB download)"
    )
    parser.add_argument("--out-dir", required=True, help="e.g. hallmark_db/")
    parser.add_argument(
        "--max-per-phrog",
        type=int,
        default=30,
        help="cap structures kept per PHROG family (default: 30, matches the validated "
        "reference build -- see docs/hallmark_db.md)",
    )
    parser.add_argument("--seed", type=int, default=0, help="stride-sample seed (unused by the stride itself, kept for reproducibility bookkeeping)")
    parser.add_argument("--foldseek", default="foldseek", help="foldseek binary/path (envs/foldseek.yaml)")
    parser.add_argument(
        "--skip-createsubdb",
        action="store_true",
        help="only write the id lists and category tables, don't invoke foldseek "
        "(useful for a quick check of the category counts before committing to the "
        "createsubdb step, which reads through the full multi-GB structure DB)",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = Path(args.phold_db_prefix)
    lookup_path = Path(f"{args.phold_db_prefix}.lookup")

    hallmark, separate = categorise(args.annots)
    print(f"hallmark PHROGs (head/connector/tail/lysis): {len(hallmark)}", file=sys.stderr)
    print(f"separate channel (integration and excision): {len(separate)}", file=sys.stderr)

    all_wanted = set(hallmark) | set(separate)
    by_phrog = load_lookup_by_phrog(lookup_path, all_wanted)

    missing = all_wanted - set(by_phrog)
    if missing:
        print(
            f"WARNING: {len(missing)} PHROG ids have no structures in the DB "
            f"(e.g. {sorted(missing)[:5]})",
            file=sys.stderr,
        )

    capped = cap_representatives(by_phrog, args.max_per_phrog, args.seed)

    def write_channel(name, phrog_ids, category_map):
        ids_path = out_dir / f"{name}_ids.tsv"
        cat_path = out_dir / f"{name}_categories.tsv"
        n_structures = 0
        with open(ids_path, "w") as ids_out, open(cat_path, "w") as cat_out:
            cat_out.write("phrog_id\tcategory\n")
            for phrog_id in sorted(phrog_ids, key=int):
                for row_idx in capped.get(phrog_id, []):
                    ids_out.write(f"{row_idx}\n")
                    n_structures += 1
                cat_out.write(f"phrog_{phrog_id}\t{category_map[phrog_id]}\n")
        print(f"{name}: {len(phrog_ids)} PHROGs, {n_structures} structures -> {ids_path}", file=sys.stderr)
        return ids_path

    # "hallmark" here matches what --hallmark-db/--hallmark-categories and
    # phables_utils/hallmark_utils.py's load_hallmark_categories() expect by
    # name (hallmark_db, hallmark_categories.tsv) -- do not rename without
    # updating those too.
    hallmark_ids_path = write_channel("hallmark", hallmark.keys(), hallmark)
    separate_ids_path = write_channel("integration_excision", separate.keys(), separate)

    if args.skip_createsubdb:
        return

    # --id-mode 0: subset by database KEY (the .lookup file's row index, column 0),
    # not by sequence content -- confirmed via `foldseek createsubdb -h`, and matches
    # the row indices load_lookup_by_phrog collected above.
    #
    # _ca (C-alpha coordinates) needs its OWN separate createsubdb call, same as
    # _ss/_h -- it is not a side effect of subsetting the main "" db. Confirmed
    # against steineggerlab/foldseek#97, which builds an AFDB subset the same
    # way: separate createsubdb calls against afdb, afdb_ss, and afdb_ca. Missed
    # in an earlier version of this script; foldseek's structural search/
    # alignment needs the coordinates, not just AA+3Di sequence, so a subDB
    # built without _ca would be incomplete.
    for name, ids_path in (("hallmark", hallmark_ids_path), ("integration_excision", separate_ids_path)):
        out_prefix = out_dir / f"{name}_db"
        for suffix in ("", "_ss", "_h", "_ca"):
            run([
                args.foldseek, "createsubdb", "--id-mode", "0",
                str(ids_path),
                f"{db_prefix}{suffix}",
                f"{out_prefix}{suffix}",
            ])

    print(
        f"\nDone. Point --hallmark-db at {out_dir / 'hallmark_db'} and "
        f"--hallmark-categories at {out_dir / 'hallmark_categories.tsv'}.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
