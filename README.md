# Multiple Sequence Alignment


## Project Description

The program computes multiple sequence alignments based on the methodology outlined in this specific research paper [Desmond G. Higgins, Paul M. Sharp, Fast and sensitive multiple sequence alignments on a microcomputer, Bioinformatics, Volume 5, Issue 2, April 1989, Pages 151–153,] (https://doi.org/10.1093/bioinformatics/5.2.151)
It leverages dynamic programming to progressively align sequences according to an initial phylogenetic tree's branching order. 

## Download this repository

```bash
git clone https://github.com/LouizaG/clustal.git
cd clustal
```
### Conda environment

Install [conda](https://docs.conda.io/en/latest/miniconda.html).

Install mamba:

```bash
conda install mamba -n base -c conda-forge
```

Create conda environment and install dependendies:

```bash
mamba env create -f binder/environment.yml
```

Load conda environment:

```bash
conda activate 3DGB
```

