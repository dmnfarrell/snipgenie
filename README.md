# Pathogenie

<img align="right" src=pathogenie/logo.png width=180px>

_Pathogenie_ is a desktop application for microbial variant calling and phylogenetic analysis from raw read data. It is primarily written to be used with bacterial isolates of MBovis but can be applied to other species. This is in early stages of development. Anyone interested in using the software is encouraged to make sugggestions on improving or adding features.

This software is written in Python and is developed with the Qt toolkit using PySide2. It was made on Ubuntu linux but is designed to also run on Windows 10 with a standalone application.

## Documentation

https://pathogenie.readthedocs.io/en/latest/

## Installation

Note for Windows users: a standalone installer will be available.

`pip install -e git+https://github.com/dmnfarrell/pathogenie.git#egg=pathogenie`

(Use pip3 on Ubuntu to ensure you use Python 3)

## Usage

Run `pathogenie` for the cli or `pathogenie-gui` for the desktop version.

## Dependencies

For Linux installs, you require Python 3 and the following packages. These will be installed automatically when using pip.

* numpy
* pandas
* matplotlib
* biopython
* pyvcf
* pyside2 (if using GUI)

Other binaries required:

* bwa
* samtools
* bcftools

The binaries can be installed with apt in Ubuntu. They are downloaded automatically in Windows.

## Screenshots

Desktop application:

<img src=img/scr1.png width=450px>

## BTBGENIE

This is software is developed as part of the **BTBGENIE** project, a DAFM funded project for development of genomic epidemiology systems for tracking and eradicating Mycobacterium bovis in Ireland.
