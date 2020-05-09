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

These binaries can be installed with apt in Ubuntu:

`sudo apt install bwa samtools bcftools`

If you want a tree to be built you should install RaXML:

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
                      symbol to split the sample names on
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
-o FILE, --outdir FILE
                      Results folder
-v, --version         Get version
-s, --test            Do test run
```

### Examples

```
snpgenie -r reference.fa -g reference.gff -i data_files -t 8 -o results
```

Add your own filters:

```
snpgenie -r reference.fa -g reference.gff -i data_files -t 8 -o results` \
 -f 'QUAL>=40 && INFO/DP>=20 && MQ>40'
```


### From Python

You can run a workflow from within Python:

```python
from sngenie import app
args = {'threads':8, 'outdir': 'results', 'labelsep':'-',
        'input':['/my/folder/',
                 '/my/other/folder'],
        'reference': None, 'overwrite':False}
W = app.WorkFlow(**args)
st = W.setup()
W.run()
```

## BTBGENIE

This is software is developed as part of the **BTBGENIE** project, a DAFM funded project for development of genomic epidemiology systems for tracking and eradicating Mycobacterium bovis in Ireland.
