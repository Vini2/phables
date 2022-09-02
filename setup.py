#!/usr/bin/env python3

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

packages = setuptools.find_packages()
package_data = {"phables_utils": ["phables_utils/*"]}

data_files = [(".", ["LICENSE", "README.md"])]

setuptools.setup(
    name="phables",
    version="0.1",
    zip_safe=True,
    author="Vijini Mallawaarachchi",
    author_email="viji.mallawaarachchi@gmail.com",
    description="PhaBles: Phage bubbles resolve bacteriophage genomes from metagenomic data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Vini2/phables",
    license="GPL-3.0",
    packages=packages,
    package_data=package_data,
    data_files=data_files,
    include_package_data=True,
    scripts=["phables"],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL-3.0 License",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "biopython",
        "python-igraph",
        "networkx",
        "scipy",
        "numpy",
        "tqdm",
        "click"
    ],
    python_requires=">=3.8",
)
