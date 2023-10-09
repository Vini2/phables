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
    python_requires=">=3.9",
    description="Phables: from fragmented assemblies to high-quality bacteriophage genomes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=get_version(),
    author="Vijini Mallawaarachchi",
    author_email="viji.mallawaarachchi@gmail.com",
    data_files=data_files,
    py_modules=["phables"],
    install_requires=[
        "snakemake>=7.14.0",
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
        "Operating System :: OS Independent",
    ],
)
