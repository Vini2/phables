import subprocess
from pathlib import Path

import pytest

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2022, Phables Project"
__license__ = "MIT"
__type__ = "Test Script"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"


TEST_ROOTDIR = Path(__file__).parent
EXEC_ROOTDIR = Path(__file__).parent.parent


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


@pytest.fixture(autouse=True)
def workingdir(tmp_dir, monkeypatch):
    """set the working directory for all tests"""
    monkeypatch.chdir(tmp_dir)


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0
    Parameters
    ----------
    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None


def test_phables(tmp_dir):
    """test phables"""
    dir_name = TEST_ROOTDIR / "data"
    graph = dir_name / "assembly_graph.gfa"
    paths = dir_name / "assembly_info.txt"
    coverage = dir_name / "edge_coverages.tsv"
    smg = dir_name / "edges.fasta.hmmout"
    phrogs = dir_name / "phrogs_annotations.tsv"

    cmd = f"python {EXEC_ROOTDIR}/phables/workflow/scripts/phables.py -g {graph} -p {paths} -b {dir_name} -hm {smg} -ph {phrogs} -c {coverage} -o {tmp_dir}"
    exec_command(cmd)


def test_combine_cov(tmp_dir):
    """test combine_cov"""
    dir_name = TEST_ROOTDIR / "script_data"
    cmd = f"python {EXEC_ROOTDIR}/phables/workflow/scripts/combine_cov.py --covpath {dir_name} --output {tmp_dir}"
    exec_command(cmd)


def test_gfa2fasta(tmp_dir):
    """test gfa2fasta"""
    dir_name = TEST_ROOTDIR / "data"
    graph = dir_name / "assembly_graph.gfa"
    cmd = f"python {EXEC_ROOTDIR}/phables/workflow/scripts/gfa2fasta.py --graph {graph} --assembler flye --output {tmp_dir}"
    exec_command(cmd)
