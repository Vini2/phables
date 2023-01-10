import os
from setuptools import setup, find_packages


def get_version():
    with open(
        os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "phables",
            "phables.VERSION",
        )
    ) as f:
        return f.readline().strip()


CLASSIFIERS = [
    "Development Status :: 3 - Alpha",
    "Programming Language :: Python :: 3",
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT license",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name="phables",
    packages=find_packages(),
    url="https://github.com/Vini2/phables",
    python_requires=">=3.8",
    description="Phables: Phage bubbles resolve bacteriophage genomes in viral metagenomic samples",
    version=get_version(),
    author="Vijini Mallawaarachchi",
    author_email="viji.mallawaarachchi@gmail.com",
    py_modules=["phables"],
    install_requires=[
        "snakemake>=7.14.0",
        "pyyaml>=6.0",
        "Click>=8.1.3",
    ],
    entry_points={
        "console_scripts": [
            "phables=phables.__main__:main"
        ]
    },
    include_package_data=True,
)
