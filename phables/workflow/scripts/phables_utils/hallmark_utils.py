from collections import defaultdict


def recover_unitig(protein_id):
    """Strip the last three underscore fields: {unitig}_{start}_{end}_{strand} -> unitig.

    Mirrors gene_utils.get_smg_unitigs exactly, so both the marker-gene and the
    hallmark-detection paths fail the same way on a naming-convention mismatch
    rather than one silently breaking while the other works.
    """
    parts = protein_id.split("_")
    return "_".join(parts[: len(parts) - 3])


def load_hallmark_categories(path):
    """phrog_id -> category string, from build_hallmark_db.py's *_categories.tsv."""
    categories = {}
    with open(path, "r") as myfile:
        next(myfile)  # header
        for line in myfile:
            phrog_id, category = line.rstrip("\n").split("\t")
            categories[phrog_id] = category
    return categories


def _parse_hits(hits_tsv):
    """Yield (unitig, phrog_id, evalue, bits) per foldseek convertalis line.

    Expects --format-output query,target,evalue,bits,fident.
    """
    with open(hits_tsv, "r") as myfile:
        for line in myfile:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            query, target, evalue, bits = fields[0], fields[1], fields[2], fields[3]
            phrog_id = target.split(":", 1)[0]  # "phrog_2:protein419484" -> "phrog_2"
            yield recover_unitig(query), phrog_id, float(evalue), float(bits)


def get_hallmark_unitigs(hits_tsv, categories, e_value, min_bits=0.0):
    """Drop-in equivalent of gene_utils.get_phrog_unitigs, foldseek-sourced.

    Returns (unitig_phrogs, phrog_dict) with the same shapes as the mmseqs path:
    unitig_phrogs[unitig] = {phrog_id, ...}; phrog_dict[phrog_id] = category string.

    component_utils.get_components does an unguarded phrog_dict[phrog] lookup for
    every id in unitig_phrogs -- safe in the mmseqs path because phrog_dict there is
    built from the entire phrog_annot.tsv, a strict superset of anything that path's
    search could return. This path's search space is narrower (foldseek can only
    return targets that exist in the hallmark DB, built from this same categories
    file), so the invariant holds by construction today -- but only if categories
    actually matches the DB that produced hits_tsv. A stale categories file from a
    previous DB build would violate that silently and surface as a KeyError deep
    inside phables, far from the real cause. Assert it here instead, at the
    boundary, where the mismatch is obvious.
    """
    unitig_phrogs = defaultdict(set)
    for unitig, phrog_id, evalue, bits in _parse_hits(hits_tsv):
        if evalue < e_value and bits >= min_bits:
            unitig_phrogs[unitig].add(phrog_id)

    hit_phrogs = {p for phrogs in unitig_phrogs.values() for p in phrogs}
    unknown = hit_phrogs - set(categories)
    if unknown:
        raise ValueError(
            f"{len(unknown)} phrog id(s) in the foldseek hits have no entry in the "
            f"categories file (e.g. {sorted(unknown)[:5]}). This means hits_tsv was "
            f"searched against a hallmark DB that doesn't match the categories file "
            f"passed in -- rebuild both from the same build_hallmark_db.py run."
        )

    return dict(unitig_phrogs), dict(categories)


def score_unitigs(hits_tsv, categories, e_value):
    """Per-unitig aggregate score, for calibration/reporting -- not required by
    get_components, which only needs presence/absence per category. "Don't just
    count hits": distinct hallmark families hit, best bitscore per family, total
    hit count. A unitig hitting three different structural families is stronger
    evidence than three hits to the same one.
    """
    per_unitig = defaultdict(lambda: defaultdict(list))  # unitig -> category -> [bits]
    for unitig, phrog_id, evalue, bits in _parse_hits(hits_tsv):
        if evalue >= e_value:
            continue
        category = categories.get(phrog_id, "unknown")
        per_unitig[unitig][category].append(bits)

    rows = []
    for unitig, by_cat in per_unitig.items():
        best_bits = {cat: max(vals) for cat, vals in by_cat.items()}
        rows.append(
            {
                "unitig": unitig,
                "n_families": len(by_cat),
                "families": ";".join(sorted(by_cat)),
                "total_hits": sum(len(v) for v in by_cat.values()),
                "best_bits_per_family": ";".join(
                    f"{cat}:{best_bits[cat]:.0f}" for cat in sorted(best_bits)
                ),
            }
        )
    return rows
