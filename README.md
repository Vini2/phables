<p align="center">
  <img src="ds-bubbles-logo.png" width="800" title="ds-bubbles logo" alt="ds-bubbles logo">
</p>

# ds-bubbles: Viral bubbles resolve bacteriophage genomes from metagenomic data

`ds-bubbles` is a tool developed to resolve bacteriophage genomes using viral bubbles in metagenomic data. `ds-bubbles` requires the the assembled contigs and assembly graph files from [Hecatomb](https://hecatomb.readthedocs.io/en/latest/). 

## Usage

You can see the following command line options of `dsbubbles.py` using `python dsbubbles.py --help`.

```
Usage: dsbubbles.py [OPTIONS]

  ds-bubbles: Resolve bacteriophage genomes from viral bubbles in metagenomic
  data.

Options:
  -g, --graph PATH          path to the assembly graph file  [required]
  -c, --contigs PATH        path to the contigs file  [required]
  -p, --paths PATH          path to the contig paths file  [required]
  -hm, --hmmout PATH        path to the contig .hmmout file  [required]
  -ph, --phrogs PATH        path to the contig phrog annotations file
                            [required]
  -ml, --minlength INTEGER  minimum length of circular contigs to consider
  -mgf, --mgfrac FLOAT      length threshold to consider single copy marker
                            genes
  -as, --alignscore FLOAT   minimum alignment score (%) for phrog annotations
  -si, --seqidentity FLOAT  minimum sequence identity for phrog annotations
  -o, --output PATH         path to the output folder  [required]
  --help                    Show this message and exit.
```
