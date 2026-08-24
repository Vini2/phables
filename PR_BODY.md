Follow-up to #68. Three independent strands: a crash fix, two flow-decomposition
speedups, and a rework of the container added in #68.

Rebases cleanly onto current `develop` (checked against `d196f5f`, including the
`QUAY_NAMESPACE` change to `build_container.yaml`).

---

## 1. Fix: phables crashes when a sample resolves no genomes

A sample where nothing is resolved and there are no unresolved phage-like edges
produces an empty `genomes_and_unresolved_edges.fasta`. CoverM 0.7.0 doesn't
print an empty table for the resulting zero-alignment BAM — it panics:

```
[WARN  coverm::contig] No primary alignments were observed for sample X
thread 'main' panicked at src/coverage_printer.rs:467:61:
index out of bounds: the len is 0 but the index is 0
```

That killed the whole run at the very last stage, after all the expensive work
had already succeeded. Hit for real on `SRR19670770`.

`coverm_bam2counts_genomes` now checks the BAM for alignments first and writes a
header-only coverage table instead of invoking CoverM when there are none.
Everything downstream already handled an empty table correctly, so the run
finishes normally with empty report tables. The header reproduces CoverM's own
exactly, including `Covered Fraction` being two words.

## 2. Performance: flow decomposition

Both changes are opt-out-safe — the second is off by default — and both were
profiled before being written rather than guessed at.

**Where the time actually goes.** For a component's MILP, *building* the model
is ~95% of the cost, not solving it (large component: 170 ms build vs 8.9 ms
solve). Two consequences, both measured: solver threads make no difference at
all (1/2/4/8 threads are flat within noise), and a solver `time_limit` does not
help either.

**a. Start the K search at a proven lower bound** (`FD_Algorithm`)

`FD_Algorithm` tried K = 1, 2, 3, … until feasible, rebuilding the whole MILP
each time. K is structural to the model, so it genuinely cannot be reused, and
flowpaths does not expose a HiGHS warm start — meaning every attempt below the
true answer was a full model build that could only return infeasible.

`get_lowerbound_k()` takes the max of the graph width and
`ceil(log2(#distinct flow values))`, both lifted from flowpaths' own
`MinFlowDecomp.get_lowerbound_k`. Both are lower bounds, so starting there
cannot skip a feasible smaller K. It costs 1–4 ms and falls back to 1 on any
error, since a lower bound is an optimisation and must never be why a component
fails to resolve.

| | |
|---|---|
| 18/18 synthetic cases | identical K, path count and path sets |
| speedup | 1.5×–4.9×, growing with component size |
| components that can't resolve within `--maxpaths` | up to 5.9× (the bound proves `K >= maxpaths` up front instead of burning the whole ladder) |

(Also annotates `data["minK"]`, which was set to a constant `2` and never read
by anything, so it isn't mistaken for the live bound.)

**b. `--mfd-workers`: run components in parallel** (default `1`, unchanged behaviour)

Components are independent, so the loop is embarrassingly parallel.
`resolve_short_parallel` chunks them, runs the **existing** `resolve_short` once
per chunk in a worker process, and merges the returned accumulators — no change
to that function's ~1400-line body, since it's already parameterised by the
component set and already returns everything it builds.

This is only sound because no component's logic depends on another's results.
That was verified against the body first: every touch of a shared accumulator is
a pure `add`/`union`/`append`, with no conditional or membership test against
them anywhere in the loop (per-component decisions use the loop-local `comp_*`
sets). It's noted in the docstring, because chunking would silently change
results if that ever stopped being true.

Chunks merge in component order, so `all_resolved_paths` is identical to the
sequential run — genomes are numbered by position, so a different order would
rename every genome without changing the biology.

| Workload | 2 | 4 | 8 workers |
|---|---|---|---|
| uniform components | 1.93× | 3.74× | 6.08× |
| realistic skew | 1.92× | 2.02× | 2.35× |

The skewed case is the one to plan around: one component was 41% of total
runtime, giving an Amdahl ceiling of 2.4× — so 2.35× is ~98% of what's
achievable. Two details mattered: more chunks than workers (so the pool can
balance an uneven workload; this alone took 8 workers from 1.92× to 2.35×), and
sending the heavy read-only inputs once per worker via a pool initializer rather
than once per chunk, so smaller chunks don't mean re-pickling the assembly graph.

Verified identical results at 2/4/8 workers, identical path ordering, and
identical output for 1 component, fewer components than workers, more components
than workers, and `workers=1`.

## 3. Container: one monolithic image, replacing the per-rule container flags

**This removes `--container` and `--prostt5-container`, both added in #68.**
Worth being explicit about, since they were merged only recently.

Those flags pointed individual rules at images via Snakemake `container:`
directives. That approach fights Snakemake: combining `container:` with `conda:`
triggers the documented "ad-hoc combination" behaviour, which builds a *fresh*
conda env *inside* the container rather than using what the image already has —
defeating the point. The rules now carry plain `conda:` directives only, and the
image satisfies them by having every per-rule env **pre-built inside it**
(hybracter's approach). Run it and there is nothing left to create at runtime,
which removes the conda-env creation race that motivated containerising at all:
each sample is its own Snakemake process, so concurrent array tasks needing the
same not-yet-built env can corrupt each other's `mamba env create`.

Also here:
- **`--gpu-backend system`** — builds no env and uses the ambient `torch` +
  `pholdlib`. Conda envs are isolated, so a rule declaring `conda:` can never see
  a torch installed outside it; this is the only way to *reuse* a known-good GPU
  torch (the container's ROCm base, or a module-loaded torch on HPC) instead of
  installing a second copy. Keeps several GB out of the image.
- `container/prebuild_envs.sh` generates its own throwaway inputs rather than
  using `tests/data/` (which is gitignored, so absent in a fresh clone or CI
  checkout), and points `phables install` at an *empty* databases dir so the
  download rules are actually in the DAG and the `curl` env gets built.
- `container/test_image.sh` fails the build if any per-rule tool is missing, if
  `torch`/`pholdlib` aren't importable, or if any pre-built env contains its own
  torch (i.e. a second copy crept in).
- `prostt5-rocm.yaml` bumped to `torch==2.9.1` (verified present on the pinned
  `rocm6.3` index, cp310–cp314).

---

## Testing

- Existing test suite passes.
- Workflow DAG verified by real `--dry-run`s across gene caller, detection mode,
  GPU backend and tree options.
- Performance and equivalence numbers above are from synthetic instances driving
  the real `FD_Inexact` / `resolve_short` code paths.
- The container built end-to-end: all 8 per-rule conda envs solve and install, no
  duplicate torch, and the pip install leaves the base image's torch untouched.

## Not verified

- No end-to-end run on a real sample through the parallel path — the equivalence
  testing used a stubbed `resolve_short` (real chunking and merging, synthetic
  per-component results).
- Speedups are from synthetic components with clean topology. Build-dominates-
  solve should hold generally since it tracks graph size, but a genuinely hard
  component would shift the ratio. The lower bound is *valid* regardless, just
  possibly looser on real graphs.
- No Setonix Apptainer run of the final image.
- Process pools copy memory per worker, so a large assembly graph at 8 workers
  may bind on RAM before CPU.

Happy to split any of the three strands into its own PR if that's easier to
review — they're independent.
