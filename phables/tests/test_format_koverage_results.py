#!/usr/bin/env python3

"""Regression test for workflow/scripts/format_koverage_results.py's parsing of
the per-genome sample_coverage.tsv.

Context
-------
postprocess.smk's per-genome coverage is produced by a direct `coverm contig`
call (rules coverm_map_genomes / coverm_bam2counts_genomes /
coverm_combine_genomes) -- koverage was removed from the workflow entirely per
PLAN.md 4.8. That table's header is:

    Sample  Contig  Count  RPKM  TPM  Mean  Covered_fraction  Variance
      0       1       2     3    4     5           6             7

It previously came from Koverage's *native* "map" mode, whose header is
different and longer:

    Sample  Contig  Count  RPM  RPKM  RPK  TPM  Mean  Median  Hitrate  Variance
      0       1       2     3    4     5    6    7      8        9        10

Under that older layout RPKM was index 4 and Mean was index 7. Getting these
two layouts mixed up silently mislabels RPKM/TPM or Mean/Variance in phables'
sample_genome_rpkm.tsv / sample_genome_mean_coverage.tsv reports rather than
erroring, so this test locks in the CoverM-mode indices the script now relies
on -- and explicitly asserts the old koverage-mode indices are NOT in use.

Run directly with:
    python3 tests/test_format_koverage_results.py
"""

import importlib.util
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT_PATH = os.path.join(
    REPO_ROOT, "workflow", "scripts", "format_koverage_results.py"
)

# Import format_koverage_results.py as a module without executing its
# `if __name__ == "__main__":` block (which requires a snakemake object).
spec = importlib.util.spec_from_file_location(
    "format_koverage_results", SCRIPT_PATH
)
format_koverage_results = importlib.util.module_from_spec(spec)
spec.loader.exec_module(format_koverage_results)


# The CoverM-mode header, as written by postprocess.smk's
# coverm_combine_genomes rule (coverm's own "<bam filename> <Method>" column
# names with the filename prefix stripped, field 0 renamed to "Contig", and a
# leading "Sample" column prepended).
COVERM_MODE_HEADER = (
    "Sample\tContig\tCount\tRPKM\tTPM\tMean\tCovered_fraction\tVariance"
)

# A corresponding data row with a distinguishable value in every column, so a
# wrong index picks up an obviously wrong number rather than a plausible one.
COVERM_MODE_ROW = (
    "sample1\tgenome_A\t100\t22.22\t44.44\t55.55\t0.9\t77.77"
)

# The older Koverage native "map" mode row this script used to receive, kept
# only so the test below can assert we are NOT still parsing with those
# indices. Verbatim shape from Koverage's coverage.smk all_sample_coverage
# rule + scripts/sampleCoverage.py.
LEGACY_KOVERAGE_MAP_MODE_ROW = (
    "sample1\tgenome_A\t100\t11.11\t22.22\t33.33\t44.44\t55.55\t66.66\t0.9\t77.77"
)


def test_parse_koverage_row_uses_coverm_mode_indices():
    strings = COVERM_MODE_ROW.split("\t")
    sample, contig, count, rpkm_val, mean_val = (
        format_koverage_results.parse_koverage_row(strings)
    )

    assert sample == "sample1", f"expected sample='sample1', got {sample!r}"
    assert contig == "genome_A", f"expected contig='genome_A', got {contig!r}"
    assert count == 100, f"expected Count=100, got {count!r}"
    # Column 3 in CoverM mode is RPKM (not TPM, which is column 4).
    assert rpkm_val == 22.22, f"expected RPKM=22.22, got {rpkm_val!r}"
    # Column 5 in CoverM mode is Mean (not Covered_fraction=6 / Variance=7).
    assert mean_val == 55.55, f"expected Mean=55.55, got {mean_val!r}"


def test_header_matches_expected_coverm_mode_layout():
    header_fields = COVERM_MODE_HEADER.split("\t")
    assert header_fields[format_koverage_results.IDX_SAMPLE] == "Sample"
    assert header_fields[format_koverage_results.IDX_CONTIG] == "Contig"
    assert header_fields[format_koverage_results.IDX_COUNT] == "Count"
    assert header_fields[format_koverage_results.IDX_RPKM] == "RPKM"
    assert header_fields[format_koverage_results.IDX_MEAN] == "Mean"


def test_legacy_koverage_indices_are_not_in_use():
    """Guard against silently reverting to the old native-"map"-mode indices.

    Parsing the CoverM-mode row with those indices would still "work" (no
    exception, both are valid floats) but would report TPM as RPKM and
    Covered_fraction as Mean -- exactly the silent mislabelling this file
    exists to prevent.
    """
    assert format_koverage_results.IDX_RPKM != 4, (
        "IDX_RPKM=4 is the legacy koverage native-map-mode index; CoverM mode "
        "puts RPKM at 3. At index 4 you would be reporting TPM as RPKM."
    )
    assert format_koverage_results.IDX_MEAN != 7, (
        "IDX_MEAN=7 is the legacy koverage native-map-mode index; CoverM mode "
        "puts Mean at 5. At index 7 you would be reporting Variance as Mean."
    )
    # And confirm the legacy row is genuinely a different shape, i.e. these two
    # layouts really are distinguishable rather than this being a no-op check.
    assert len(LEGACY_KOVERAGE_MAP_MODE_ROW.split("\t")) == 11
    assert len(COVERM_MODE_ROW.split("\t")) == 8


def main():
    tests = [
        test_parse_koverage_row_uses_coverm_mode_indices,
        test_header_matches_expected_coverm_mode_layout,
        test_legacy_koverage_indices_are_not_in_use,
    ]
    failures = 0
    for test in tests:
        try:
            test()
            print(f"PASS: {test.__name__}")
        except AssertionError as e:
            failures += 1
            print(f"FAIL: {test.__name__}: {e}")

    if failures:
        print(f"\n{failures} test(s) failed")
        sys.exit(1)
    else:
        print("\nAll tests passed")


if __name__ == "__main__":
    main()
