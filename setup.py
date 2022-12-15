#!/usr/bin/env python3

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

packages = setuptools.find_packages()
package_data = {"phables_utils": ["phables_utils/*", "phables_utils/support/*", "phables_utils/phrogs/*"]}

data_files = [(".", ["LICENSE", "README.md"])]

setuptools.setup(
    name="phables",
    version="0.1.0a3",
    zip_safe=True,
    author="Vijini Mallawaarachchi",
    author_email="viji.mallawaarachchi@gmail.com",
    description="Phables: Phage bubbles resolve bacteriophage genomes in viral metagenomic samples",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Vini2/phables",
    license="MIT",
    packages=packages,
    package_data=package_data,
    data_files=data_files,
    include_package_data=True,
    scripts=["phables"],
    entry_points = {
        'console_scripts': [
            'combine_cov=phables_utils.support.combine_cov:main',
            'gfa2fasta=phables_utils.support.gfa2fasta:main'
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "biopython",
        "python-igraph",
        "networkx>=2.8.6",
        "scipy",
        "numpy",
        "pandas",
        "tqdm",
        "click",
        "pysam",
        "more-itertools",
        "gurobipy",
    ],
    python_requires=">=3.8",
)
