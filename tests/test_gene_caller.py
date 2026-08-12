import sys
from pathlib import Path

import pytest

TEST_ROOTDIR = Path(__file__).parent
EXEC_ROOTDIR = Path(__file__).parent.parent
SCRIPTS_DIR = EXEC_ROOTDIR / "phables" / "workflow" / "scripts"
TEST_DATA = EXEC_ROOTDIR / "phables" / "test_data"

sys.path.insert(0, str(SCRIPTS_DIR))

from gene_caller import call_genes, gene_id, read_fasta, write_proteins  # noqa: E402

sys.path.insert(0, str(SCRIPTS_DIR / "phables_utils"))


def strip_three(name):
    """Mirror of gene_utils.get_smg_unitigs unitig-name recovery."""
    parts = name.split("_")
    return "_".join(parts[: len(parts) - 3])


@pytest.fixture(scope="module")
def edges_fasta(tmp_path_factory):
    """Build a unitig FASTA from the bundled test assembly graph."""
    path = tmp_path_factory.mktemp("data") / "edges.fasta"
    records = []
    with open(TEST_DATA / "assembly_graph.gfa") as handle:
        for line in handle:
            if line.startswith("S"):
                fields = line.rstrip("\n").split("\t")
                records.append(f">{fields[1]}\n{fields[2]}")
    path.write_text("\n".join(records) + "\n")
    return path


@pytest.fixture(scope="module")
def called(edges_fasta):
    records = list(read_fasta(edges_fasta))
    return records, call_genes(records, threads=1)


def test_reads_all_segments(edges_fasta, called):
    records, _ = called
    assert len(records) == 16
    assert records[0][0] == "edge_1"


def test_genes_are_called(called):
    _, results = called
    assert sum(len(genes) for _, genes in results) > 0


def test_header_is_fraggenescan_compatible(called):
    """Ids must be {unitig}_{start}_{end}_{strand}, else get_smg_unitigs breaks."""
    _, results = called
    for unitig, genes in results:
        for gene in genes:
            ident = gene_id(unitig, gene)
            assert ident.startswith(f"{unitig}_")
            suffix = ident[len(unitig) + 1 :].split("_")
            assert len(suffix) == 3
            assert suffix[0].isdigit() and suffix[1].isdigit()
            assert suffix[2] in ("+", "-")


def test_unitig_name_survives_round_trip(called):
    """The regression that motivated the header choice.

    Pyrodigal's default `{seqid}_{index}` naming makes strip_three return "" for
    every protein, which silently empties the marker-gene set and disables
    bacterial filtering entirely.
    """
    _, results = called
    for unitig, genes in results:
        for gene in genes:
            assert strip_three(gene_id(unitig, gene)) == unitig


def test_default_naming_would_regress(called):
    """Guard the above: prove the naive naming really does collapse to ''."""
    _, results = called
    unitig, genes = next((u, g) for u, g in results if len(g) > 0)
    assert strip_three(f"{unitig}_1") == ""


def test_partial_genes_are_kept(called):
    """Unitigs are graph fragments; FragGeneScan ran with -complete=0 for this."""
    _, results = called
    partial = sum(
        gene.partial_begin + gene.partial_end for _, genes in results for gene in genes
    )
    assert partial > 0


def test_translations_have_no_trailing_stop(called, tmp_path):
    """FragGeneScan does not emit stop codons; HMMER input should match."""
    _, results = called
    out = tmp_path / "proteins.faa"
    write_proteins(results, out)
    for line in out.read_text().splitlines():
        if not line.startswith(">"):
            assert not line.endswith("*")


def test_output_is_deterministic(called, tmp_path):
    _, results = called
    first, second = tmp_path / "a.faa", tmp_path / "b.faa"
    write_proteins(results, first)
    write_proteins(results, second)
    assert first.read_text() == second.read_text()
