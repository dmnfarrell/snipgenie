Introduction
============

**snpgenie** is a desktop and command line tool for microbial variant calling and phylogenetic analysis from raw read data. It is primarily written to be used with bacterial isolates of MBovis but can be applied to other species. This is in early stages of development. Anyone interested in using the software is encouraged to make sugggestions on improving or adding features.

This software is written in Python and is developed with the Qt toolkit using PySide2. It was made on Ubuntu linux but is designed to also run on Windows 10 with a standalone application.

Command line tool
-----------------

This tool works from the command line and via Python scripts. Unlike many other SNP calling pipelines, it is also designed to have a graphical user interface, which is in development.

Current Features
----------------

* load multiple fastq files and process together
* view fastq quality statistics
* trim reads
* align to reference
* view bam alignments
* call variants
* filter variants
* create SNP core multiple sequence alignment
* create phylogenetic tree

Links
-----

http://dmnfarrell.github.io/snpgenie

Installation
============

Linux
-----

With pip::

  pip install -e git+https://github.com/dmnfarrell/snpgenie.git#egg=snpgenie

Install binary dependencies::

  sudo apt install bcftools samtools bwa

Windows
-------

A standalone installer will be used to deploy on windows.

Mac OSX
-------

Not tested. You can try the Linux instructions possibly with bioconda for the binaries.
