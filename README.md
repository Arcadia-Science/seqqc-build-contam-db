# Build a sourmash database to detect common sources of contamination in new sequencing data

**Table of Contents**
* [Introduction](#introduction)
* [Accessing the built database](#accessing-the-built-database)
* [Adding new data to the database](#adding-new-data-to-the-database)
* [Getting started with this repository](#getting-started-with-this-repository)
* [Running this repository on AWS](#running-this-repository-on-AWS)

## Introduction

This repository builds a database that can be used with `sourmash gather` to detect contamination in FASTQ files.
The database is used as part of the [seqqc pipeline](https://github.com/Arcadia-Science/seqqc) for rapid quality control of new sequencing data.
See the [seqqc pub here](https://arcadia-research.pubpub.org/pub/resource-seqqc) (DOI: 10.57844/arcadia-cxn6-ch62).

Contamination in sequencing data can come from a lot of different sources.
Below we outline six common types of contamination. 
For the first five, we also outline the identification approach incorporated into the database built by this repository.

1. **Contam type 1**: contamination from barcode/index hopping. This happens most frequently for low-biomass samples and is an illumina artifact. 
     **Computational identification approach**: Screen for sequences of model organisms that are sequenced frequently, as these will be the sequences that are most likely to occur as contaminants because they are the most likely things to be sequenced at any given time. If we get mouse in our metagenome or in chlammy rna seq (esp before we have a terrarium), it’s probably from barcode hopping
2. **Contam type 2**: contamination from humans handling the sample. This could be human sequence, or sequence from microbes that live on human skin/oral cavity (like S. aureus).
    **Computational identification approach**: Include human DNA and human skin/oral microbiome species in the database.
4. **Contam type 3**: kit contamination. Kits and reagents have their own microbiome and so DNA extracted from these organisms can sneak into the sample
    **Computational identification approach**: Add most common kit contaminant organisms to the database. This paper reviews and lists common kit contaminants: 
    > Eisenhofer, R., Minich, J. J., Marotz, C., Cooper, A., Knight, R., & Weyrich, L. S. (2018). Contamination in Low Microbial Biomass Microbiome Studies: Issues and Recommendations. Trends in Microbiology. doi:10.1016/j.tim.2018.11.003
6. **Contam type 4**: there’s lab contamination, so accidentally extracting DNA or RNA from other organisms that are in the lab. This would be chlammy for Arcadia, and any other organism that is brought into the lab
    **Computational identification approach**: Select species Arcadian's work with, create signatures for those (masked) genomes, and add them to the database
5. **Contam type 5**: spike in contamination. Illumina spikes phiX into many of of its sequencing runs. 
    **Computational identification approach**: create a sourmash signature for phix and include it in the contamination data base
8. **Contam type 6**: Contamination in the sample material itself. This might be something like mycorrhizal fungi that sneaks into a plant genome sequencing run.
    **Computational identification approach**: The best way to identify this type of contamination would be to screen with all known genomic/transcriptomic data, but this doesn't scale super well. This pipeline might do a poor job of detecting this type of contamination.

## Accessing the built database

The database is available in [this OSF repository](https://osf.io/sndz5/) (DOI: 10.17605/OSF.IO/SNDZ5).

## Adding new data to the database 

Two files are used to determine which genomes are included in the final database:

* `inputs/genome_contaminants.csv`: Genome accessions for genomes to include in the database. Accessions are GenBank or RefSeq format (`GCF_*` or `GCA_*`). This Snakefile uses the `CDS_from_genomic.fna` files, so if a genome is missing this file, it cannot be included.
* `inputs/doi10.1016j.tim.2018.11.003-table1.csv`: Genuses of bacteria or archaea that should be included in the database. To keep the final database small, one genome from each species under the specified genus is extracted from the GTDB database.

The files in this repository reflect the genomes that were included in the database posted on OSF.
To include new genomes, add to these files and re-run the snakefile (see instructions below).

## Getting started with this repository

This repository uses snakemake to run the pipeline and conda to manage software environments and installations.
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).
After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.
We installaed Miniconda3 version `py39_4.12.0` and mamba version `0.15.3`.

```
mamba env create -n seqqcdb --file environment.yml
conda activate seqqcdb
```

To start the pipeline, run:
```
snakemake --use-conda -j 2
```

### Running this repository on AWS

This repository was executed on an AWS EC2 instance (Ubuntu 22.04 LTS ami-085284d24fe829cd0, t2.2xlarge, 200 GiB EBS `gp2` root storage).
The instance was configured using the following commands:

```
curl -JLO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh # download the miniconda installation script
bash Miniconda3-latest-Linux-x86_64.sh # run the miniconda installation script. Accept the license and follow the defaults.
source ~/.bashrc # source the .bashrc for miniconda to be available in the environment

# configure miniconda channel order
conda config --add channels defaults 
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install mamba # install mamba for faster software installation.
```
