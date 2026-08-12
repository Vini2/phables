import os

from setuptools import find_packages, setup


def get_version():
    with open(
        os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "phables",
            "phables.VERSION",
        )
    ) as f:
        return f.readline().strip()


with open("README.md", "r") as fh:
    long_description = fh.read()


data_files = [(".", ["LICENSE", "README.md"])]

setup(
    name="phables",
    packages=find_packages(),
    url="https://github.com/Vini2/phables",
    # >=3.11 because snakemake>=8 (below) requires it -- every 8.x/9.x release on
    # PyPI declares python_requires>=3.11, confirmed via `pip install snakemake>=8`
    # against a 3.10 env, which refused to resolve any 8.x/9.x version at all.
    python_requires=">=3.11",
    description="Phables: from fragmented assemblies to high-quality bacteriophage genomes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=get_version(),
    author="Vijini Mallawaarachchi",
    author_email="viji.mallawaarachchi@gmail.com",
    data_files=data_files,
    py_modules=["phables"],
    install_requires=[
        # >=8: Snakemake 8.0.0 introduced the external executor-plugin architecture
        # (snakemake-executor-plugin-slurm etc) and the set-resources/executor: keys
        # in workflow profiles -- 7.x's older built-in --cluster/--drmaa flags don't
        # support per-rule resource auto-translation (e.g. resources: gpu=1 ->
        # --gpus=1) at all. This is what predict_3di/scan_hallmark's optional GPU
        # routing (phables/workflow/profiles/gpu-routing/) needs. Tested against
        # 9.24.0 specifically (the current latest at the time this was written) --
        # not against an early 8.0.x release, which Snakemake's own changelog
        # cautions may have rough edges ("you might initially experience bugs and
        # glitches"), so lean towards a recent 8.x/9.x rather than the floor.
        "snakemake>=8.0",
        "pyyaml>=6.0",
        "click>=8.1.3",
        "metasnek>=0.0.5",
        "snaketool-utils>=0.0.4",
    ],
    entry_points={"console_scripts": ["phables=phables.__main__:main"]},
    include_package_data=True,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Operating System :: MacOS",
        "Operating System :: POSIX",
    ],
)
