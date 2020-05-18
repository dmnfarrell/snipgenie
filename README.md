[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Documentation Status](https://readthedocs.org/projects/snpgenie/badge/?version=latest)](https://snpgenie.readthedocs.io/en/latest/?badge=latest)

# snpgenie

<img align="right" src=snpgenie/logo.png width=180px>

_snpgenie_ is a tool for microbial variant calling and phylogenetic analysis from raw read data. It was primarily written to be used with bacterial isolates of M. Bovis but can be applied to other species. You need a good quality reference genome to align to. Anyone interested in using the software is encouraged to make suggestions on improving or adding features.

This software is written in Python. It was developed on Ubuntu linux but is designed to also run on Windows 10 with a standalone application. The GUI is made using the Qt toolkit using PySide2.

## Current Features

* load multiple fastq files and process together
* view fastq quality statistics
* trim reads
* align to reference
* view bam alignments
* call variants
* filter variants
* create SNP core multiple sequence alignment
* create phylogenetic tree

## Installation

Note for Windows users: a GUI with standalone installer will be available shortly.

`pip install -e git+https://github.com/dmnfarrell/snpgenie.git#egg=snpgenie`

Notes: You may need to use pip3 on Ubuntu to ensure you use Python 3. Use sudo if installing system-wide.

## Dependencies

For Linux installs, you require Python 3 and the following packages. These will be installed automatically when using pip.

* numpy
* pandas
* matplotlib
* biopython
* pyvcf
* pyfaidx
* pyside2 (if using GUI)

Other binaries required:

* bwa
* samtools
* bcftools
* tabix
* parallel

These binaries can be installed with apt in Ubuntu:

`sudo apt install bwa samtools bcftools tabix parallel`

If you want a tree to be built you should install RaXML, but it's optional:

`sudo apt install raxml`

The binaries are downloaded automatically in Windows.

## Usage

Run `snpgenie` for the cli or `snpgenie-gui` for the desktop version. You require a reference genome and reads in fastq format at minimum as input.

### Command line options

This will run the entire process based on a set of options given at the terminal::
```
-h, --help            show this help message and exit
-i FILE, --input FILE
                      input folder(s)
-e LABELSEP, --labelsep LABELSEP
                      symbol to split the sample labels on
-r FILE, --reference FILE
                      reference genome filename
-g FILE, --gff FILE   reference gff, optional
-w, --overwrite       overwrite intermediate files
-m, --trim            whether to trim fastq files
-q QUALITY, --quality QUALITY
                      trim quality
-f FILTERS, --filters FILTERS
                      variant calling post-filters
-t THREADS, --threads THREADS
                      cpu threads to use
-b, --buildtree       whether to try to build a phylogenetic tree
-o FILE, --outdir FILE
                      Results folder
-v, --version         Get version
-d, --dummy           Setup samples but don't run
```

### Examples

```
snpgenie -r reference.fa -g reference.gff -i data_files -o results
```

Add your own filters and provide threads:

```
snpgenie -r reference.fa -g reference.gff -i data_files -t 8 -o results` \
 -f 'QUAL>=40 && INFO/DP>=20 && MQ>40'
```

## Inputs

Folders are searched recursively for inputs with extensions `*.f*q.gz`. So be careful you don't have files in the folders you don't want included. So the following file structure will load both sets of files if you provide the parent folder as input.

```
data/
├── ERR1588781
│   ├── ERR1588781_1.fq.gz
│   └── ERR1588781_2.fq.gz
└── ERR1588785
    ├── ERR1588785_1.fastq.gz
    └── ERR1588785_2.fastq.gz
```

Filenames are parsed and a sample name extracted for each pair (if paired end). This is simply done by splitting on the _ symbol. So a file called /path/13-11594_S85_L001-4_R1_001.fastq.gz will be given a sample name 13-11594. As long as the sample names are unique this is ok. If you had a file names like A_2_L001-4_R1_001,  A_3_L001-4_R1_001 you should split on '-' instead. You can specify this in the labelsep option.

## Use from Python

You can run a workflow from within Python by importing the snpgenie package and invoking the `WorkFlow` class. You need to provide the options in a dictionary with the same keywords as the command line. Notice in this example we are loading files from two folders.

## FAQ

_The run was stopped during execution, can it be resumed?_

Yes, by default the program won't overwrite intermediate files when re-run. So just run it again. Make sure there are no old tmp.****.bam files in the mapped folder if an alignment got interrupted.

_My sample files are not being parsed properly._

This may be because your sample names are unusual. The program extracts the unique sample names from the files by using the '_' symbol as delimeter. If your names differ you can supply a different delimeter with the labelsep option.

## BTBGENIE

This is software is developed as part of the **BTBGENIE** project, a DAFM funded project for development of genomic epidemiology systems for tracking and eradicating Mycobacterium bovis in Ireland.
