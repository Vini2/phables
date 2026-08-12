# The hallmark structure database

`--phagedetection prostt5-foldseek` needs two inputs: a Foldseek structure
subDB of phage hallmark proteins (`--hallmark-db`) and its matching
PHROG-category table (`--hallmark-categories`). `phables install` fetches a
pre-built copy of both automatically, the same way it already fetches
`marker.hmm` and the PHROGs MMseqs profile DB — most users never need
anything on this page beyond knowing it happens. It exists for anyone who
wants to understand what's actually in that download, verify it, or rebuild
it themselves (a newer phold structure DB release, a different category
selection, etc.).

## What "hallmark" means here

`phables_utils/component_utils.py`'s `get_components` screens phage-like
components for structural/virion evidence — specifically the PHROG functional
categories **head and packaging**, **connector**, **tail**, and **lysis**.
Those four categories are what "hallmark" refers to throughout this feature;
everything else in the PHROG catalogue (metabolism, DNA/RNA processing,
unknown function, ...) is irrelevant to this specific structural-evidence
check and is dropped.

**Integration and excision** (integrases, recombinases) is kept as a
*separate* id list/subDB rather than folded into the hallmark set.
Integrase-family structures are a real lysogeny signal, but integrases and
recombinases are abundant on bacterial chromosomes and mobile genetic elements
too — mixing them into the structural-hallmark evidence would make that
signal noisier, not stronger.

## Source data: phold's structure database

The subDB is built by subsetting the **full phold structure database**
([`gbouras13/phold`](https://github.com/gbouras13/phold)) — not by predicting
structures from scratch. Two files from that download are needed:

- `all_phold_structures` (the Foldseek structure DB itself, plus its
  `_ss`/`_h` companion files and `.lookup` index — ~1.36M entries, several GB)
- `phold_annots.tsv` (the PHROG annotation table mapping each structure to a
  PHROG id, product name, and functional category)

`phold_annots.tsv` uses **CRLF line endings** — every field comparison in
`build_hallmark_db.py` strips a trailing `\r` explicitly. Skipping this
silently matches zero rows instead of erroring, which is exactly the failure
mode that was hit building the reference DB below — worth knowing if you ever
modify this script.

## Why a representative cap, and why a stride sample

Some PHROG families have far more structures in phold's DB than others.
Keeping every single one would bloat the subDB for no real benefit (Foldseek
search quality saturates well before "every known example"), so
`--max-per-phrog` (default 30) caps each PHROG family's contribution.

The cap is applied via a **deterministic stride sample** over the family's
sorted `.lookup` row indices — not the first N. Row order in phold's `.lookup`
file tends to cluster by source dataset, so taking the first N would bias
toward whatever dataset happened to be indexed first for that PHROG, rather
than spreading the cap across the diversity that actually exists in the
family. A stride (`ordered[int(i * len(ordered)/max_per_phrog)]` for
`i in range(max_per_phrog)`) picks evenly across the whole sorted range
instead.

## Rebuilding it

Not a normal-use step — `phables install` already fetches a pre-built copy
(see below). Rebuild only if you want a newer phold structure DB snapshot, a
different `--max-per-phrog` cap, or a different category selection.

```bash
mamba env create -f phables/workflow/envs/foldseek.yaml -n foldseek
conda activate foldseek

python phables/workflow/scripts/build_hallmark_db.py \
    --phold-db-prefix /path/to/all_phold_structures \
    --annots /path/to/phold_annots.tsv \
    --out-dir hallmark_db/ \
    --max-per-phrog 30
```

This writes, under `hallmark_db/`:

| File | What it is |
|---|---|
| `hallmark_db*` (`hallmark_db`, `hallmark_db_ss`, `hallmark_db_h`, `hallmark_db_ca` + `.index`/`.dbtype`) | The Foldseek subDB — pass its prefix (`hallmark_db/hallmark_db`) to `--hallmark-db`. `_ca` (C-alpha coordinates) needs its own `createsubdb` call same as `_ss`/`_h` — not a side effect of the main `""` call, confirmed against [steineggerlab/foldseek#97](https://github.com/steineggerlab/foldseek/issues/97) |
| `hallmark_categories.tsv` | `phrog_id\tcategory` — pass to `--hallmark-categories` |
| `hallmark_ids.tsv` | The `.lookup` row indices `createsubdb` was given (intermediate; not needed at run time) |
| `integration_excision_db*`, `integration_excision_categories.tsv`, `integration_excision_ids.tsv` | The separate integration/excision channel discussed above — not currently wired into `--hallmark-db` (phables only consumes the hallmark channel today), kept in case a future check wants it |

Pass `--skip-createsubdb` first if you just want to sanity-check the category
counts before committing to the `createsubdb` step, which reads through the
full multi-GB structure DB and is the slow part.

**`createsubdb` leaves `.lookup`/`.source` as symlinks back to the *original*
phold DB's own files** (not copies) — expected: `--id-mode 0` subsets by the
original DB's own row keys rather than renumbering, so those original files
remain valid for resolving them, and foldseek doesn't bother duplicating
potentially-huge files it doesn't need to. **Confirmed by real testing that
phables doesn't need either file at query time** — `foldseek convertalis`
resolves the `target`/`query` columns via the `_h` (header) database, which
`createsubdb` builds properly for the subset (not symlinked); `.lookup`/
`.source` are only consulted by `createsubdb` itself, as a *build-time* input,
not by anything downstream that reads the resulting subDB. Safe (and
recommended) to leave both out of a packaged tarball — see below.

## Reference build — real numbers, so you can sanity-check your own

Built once against a real, complete phold structure DB download (not a
subsample) as part of validating this feature end-to-end:

- **5,389** hallmark PHROGs (head and packaging + connector + tail + lysis)
- **436** integration-and-excision PHROGs
- **49,869** hallmark structures after capping (`--max-per-phrog 30`)
- **2,610** integration-and-excision structures after capping
- **~138 MB** total subDB size (real, non-symlinked file content) —
  comfortably inside "low hundreds of thousands of structures,
  page-cacheable" territory
- Verified **queryable**, not just built: a self-search of the hallmark subDB
  returned biologically sensible cross-hits (`phrog_2` ↔ `phrog_5653`, both
  independently annotated "terminase large subunit / head and packaging" in
  `phold_annots.tsv`)

If your own build's PHROG-category counts differ substantially from the first
two numbers above, that's worth investigating before trusting the result —
those two counts depend only on `phold_annots.tsv`'s categories, not on which
structures happen to be in your particular phold DB snapshot, so they should
be stable across phold DB versions.

## Packaging a rebuild for redistribution

If you rebuild and want `phables install` to fetch your new version (updating
`hallmark_db_url` in `phables/config/databases.yaml`), package it from
*inside* the output directory, excluding `.lookup`/`.source` (build-time-only,
confirmed unneeded at query time above) and — check first — anything left
over from a previous packaging attempt in the same directory:

```bash
cd hallmark_db/
ls   # check there's nothing stray in here before archiving everything with `.`
tar --exclude='*.lookup' --exclude='*.source' -czf hallmark_db.tar.gz .
```

The real, currently-hosted build (`hallmark_db_url` below) skipped the
`--exclude` and instead just never had `.lookup`/`.source` present when it was
packaged — either is fine, the goal is simply that neither ends up in the
tarball. Worth knowing either way, since a `tar czf ... .` run from inside a
directory that already has a `.tar.gz` sitting in it from a previous attempt
will happily include that stale archive inside the new one — harmless at
runtime (nothing reads a `.tar.gz` member), but worth a clean `ls` first.

`install.smk`'s `hallmark_db_download` rule extracts this tarball's members
at its own root into a fresh `hallmark_db/` directory it creates itself
(`mkdir -p {output}; tar -xf {file} -C {output}`) — matching exactly the `cd
hallmark_db/ && tar ... .` layout above. A tarball built any other way (e.g.
with a top-level `hallmark_db/` folder already inside it) will end up
double-nested on extraction.

## Why phold's full structure DB isn't bundled

`phables install` fetches the *pre-built subDB* (~105MB, hosted on
[Zenodo](https://zenodo.org/records/21884331)) automatically, the same way
it fetches `marker.hmm`/the PHROGs MMseqs profile DB — but it doesn't fetch
or depend on phold's own full structure database (the multi-GB,
separately-versioned source this subDB was built from) at install time.
Rebuilding only matters if you want a newer phold snapshot or different
build parameters than the shipped reference build above; get phold's
structure DB yourself for that (see phold's own docs for current download
instructions).
