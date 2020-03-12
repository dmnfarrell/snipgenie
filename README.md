# btbgenie

The **BTBGENIE** project is a DAFM funded project for development of genomic epidemiology systems for tracking and eradicating Mycobacterium bovis in Ireland.

## Pathogenie

<img align="right" src=pathogenie/pathogenie/logo.png width=180px>

_Pathogenie_ is a desktop application for microbial variant calling and phylogenetic analysis from raw read data. It is primarily written to be used with bacterial isolates of MBovis but can be applied to other species. This is in early stages of development. Anyone interested in using the software is encouraged to make sugggestions on improving or adding features.

This software is written in Python and is developed with the Qt toolkit using PySide2. It was made on Ubuntu linux but is designed to also run on Windows 10 with a standalone application.

## Installation

```
git clone https://github.com/dmnfarrell/btbgenie.git
```

Then move to the directory and run
```
python -m pathogenie.gui
```

Will be made available on pip at a future date.

## Dependencies

You require Python 3 and the following packages:

* numpy
* pandas
* matplotlib
* biopython

Other binaries required:

* bwa
* bamtools

This can be installed with apt in Ubuntu.

## Screenshots

<img src=img/scr1.png width=450px>
