Introduction
============

**snipgenie** is a desktop and command line tool for microbial variant calling and phylogenetic analysis from raw read data. It is primarily written to be used with bacterial isolates of MBovis but can be applied to other species. This is in early stages of development. Anyone interested in using the software is encouraged to make sugggestions on improving or adding features.

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

* https://github.com/dmnfarrell/snipgenie

Installation
============

Linux
-----

With pip::

  pip install -e git+https://github.com/dmnfarrell/snipgenie.git#egg=snipgenie

Note: You may need to use pip3 on Ubuntu to ensure you use Python 3.
Use sudo if installing system-wide. Running this also requires you have git installed.
The same command can be used to update to the latest version.

Install binary dependencies::

  sudo apt install bcftools samtools bwa

Windows
-------

The pip instructions will work if you have installed Python for Windows.
A standalone installer will be used to deploy on windows.

Dependencies
============

For Linux installs, you require Python 3 and the following packages.
These will be installed automatically when using pip.
::

  numpy
  pandas
  matplotlib
  biopython
  pyvcf
  pyfaidx
  pyside2 (GUI only)
  toytree (GUI only)

Other binaries required::

  bwa
  samtools
  bcftools
  tabix
  parallel

These binaries can be installed with apt in Ubuntu::

  sudo apt install bwa samtools bcftools tabix parallel

If you want a tree to be built you should install RaXML, but it's optional::

  sudo apt install raxml

The binaries are downloaded automatically in Windows.

Funding
=======

The development of this software was largely enabled through funding by the Irish Department of Agriculture.
