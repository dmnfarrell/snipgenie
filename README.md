[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Documentation Status](https://readthedocs.org/projects/snipgenie/badge/?version=latest)](https://snipgenie.readthedocs.io/en/latest/?badge=latest)

# SNiPgenie

<img align="right" src=snipgenie/logo.svg width=180px>

_SNiPgenie_ is a tool for microbial variant calling and phylogenetic analysis from raw read data. It was primarily written to be used with bacterial isolates of M. bovis but can be applied to other species. You need a good quality reference genome to align to. Anyone interested in using the software is encouraged to make suggestions on improving or adding features.

This software is written in Python. It was developed on Ubuntu linux but can also run on Windows 10 via WSL. The GUI is made using the Qt toolkit using PyQt5/PySide2.

## Current Features

* load multiple fastq files and process together
* view fastq quality statistics
* trim reads
* align to reference (fasta file with single chromosome)
* view bam alignments
* call variants
* filter variants
* create SNP core multiple sequence alignment
* create phylogenetic tree

## Installation

### Linux

`pip install -e git+https://github.com/dmnfarrell/snipgenie.git#egg=snipgenie`

Notes: You may need to use pip3 on Ubuntu to ensure you use Python 3. Running this also requires you have **git** installed. The same command can be used to update to the latest version.

### WSL on Windows

The tool can be run on Windows using Windows Subsystem for Linux (WSL). This is [easy to install](https://www.omgubuntu.co.uk/how-to-install-wsl2-on-windows-10) from the Micrsoft store. (You can also run `wsl --install' from a terminal window). You then install Ubuntu 22.04 LTS from the store. Just make sure your version of Windows 10 is also reasonably up to date.

Open Ubuntu from the start menu and run it. It will bring up a terminal. This is just like being in normal Linux. You should first install some required packages:

```sudo apt install python3-pip x11-apps libqt5extras5```

You can then install the app using pip as above.

## Dependencies

For Linux installs, you require Python 3 and the following packages. These will be installed automatically when using pip.

* numpy
* pandas
* matplotlib
* biopython
* pyvcf
* pyfaidx
* pyqt5 (GUI only)
* toytree (GUI only)

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

Run `snipgenie` for the cli or `snipgenie-gui` for the desktop version. You require a reference genome and reads in fastq format at minimum as input.

### Command line options

This will run the entire process based on a set of options given at the terminal::
```
  -h, --help            show this help message and exit
  -i FILE, --input FILE
                        input folder(s)
  -M FILE, --manifest FILE
                        manifest file with samples, optional
                        overrides input
  -e LABELSEP, --labelsep LABELSEP
                        symbol to split the sample labels on
  -x LABELINDEX, --labelindex LABELINDEX
                        position to extract label in split filenames
  -r FILE, --reference FILE
                        reference genome filename
  -S SPECIES, --species SPECIES
                        set the species reference genome, overrides -r.
                        possible values are
                        Mbovis-AF212297, MTB-H37Rv, MAP-K10,
                        M.smegmatis-MC2155,
                        Mycoplasmabovis-PG45, Sars-Cov-2
  -g FILE, --genbank_file FILE
                        annotation file, optional
  -t THREADS, --threads THREADS
                        cpu threads to use
  -w, --overwrite       overwrite intermediate files
  -T, --trim            whether to trim fastq files
  -U, --unmapped        whether to save unmapped reads
  -Q QUALITY, --quality QUALITY
                        right trim quality, default 25
  -f FILTERS, --filters FILTERS
                        variant calling post-filters
  -m MASK, --mask MASK  mask regions from a bed file
  -c, --custom          apply custom filters
  -p PLATFORM, --platform PLATFORM
                        sequencing platform, change to ont if
                        using oxford nanopore
  -a ALIGNER, --aligner ALIGNER
                        aligner to use, bwa, subread, bowtie or minimap2
  -b, --buildtree       whether to build a phylogenetic tree, requires RaXML
  -N BOOTSTRAPS, --bootstraps BOOTSTRAPS
                        number of bootstraps to build tree
  -o FILE, --outdir FILE
                        Results folder
  -q, --qc              Get version
  -d, --dummy           Check samples but don't run
  -X, --test            Test run
  -v, --version         Get version
```

### Examples

Call with your own reference fasta file:

```
snipgenie -r reference.fa -i data_files -o results
```

Use an in built species genome as reference. This will also supply an annotation file. The current options are `Mbovis-AF212297, MTB-H37Rv, MAP-K10, M.smegmatis-MC2155`:

```
snipgenie -S Mbovis-AF212297 -i data_files -o results
```

Provide more than one folder:

```
snipgenie -r reference.fa -i data_files1 -i data_files2 -o results
```

Provide an annotation (genbank format) for consequence calling:

```
snipgenie -r reference.fa -g reference.gb -i data_files -o results
```

Add your own filters and provide threads:

```
snipgenie -r reference.fa -i data_files -t 8 -o results` \
 -f 'QUAL>=40 && INFO/DP>=20 && MQ>40'
```

### Aligners

You can use any one of the following aligners: bwa, subread, bowtie or minimap2. These should be present on your system, unless using the Windows version. Note that for oxford nanopore reads you should use minimap2 and specify the platform as 'ont'.

### Mask file

You can selectively mask snp sites such as those contained in transposons or repetitive regions from being included in the output. You need to provide a bed file with the following columns: chromosome name, start and end coordinates of the regions. There is currently a built-in mask file used for M.bovis and of you select this genome as reference using the --species option it will be used automatically.

```
LT708304.1 	 105359 	 106751
LT708304.1 	 131419 	 132910
LT708304.1 	 149570 	 151187
LT708304.1 	 306201 	 307872
```

## Inputs

You can provide single folder with all the files in one place or multiple folders. Folders are searched recursively for inputs with extensions `*.f*q.gz`. So be careful you don't have files in the folders you don't want included. The following file structure will load both sets of files if you provide the parent folder as input. You can also provide multiple separate folders using -i as shown above.

For example if you provide -i data with the following structure:

```
data/
├── ERR1588781
│   ├── ERR1588781_1.fq.gz
│   └── ERR1588781_2.fq.gz
└── ERR1588785
    ├── ERR1588785_1.fastq.gz
    └── ERR1588785_2.fastq.gz
```

Filenames are parsed and a sample name is extracted for each pair (if paired end). This is simply done by splitting on the _ symbol by default. So a file called /path/13-11594_S85_L001-4_R1_001.fastq.gz will be given a sample name 13-11594. As long as the sample names are unique this is ok. If you had a file names like A_2_L001-4_R1_001, A_3_L001-4_R1_001 you should split on '-' instead. You can specify this in the `labelsep` option. The workflow won't run unless sample names are unique.

### Manifest file

You can use a manifest file (-M) with the sample names and files in a table if parsing folders won't work for you. This could be useful if your files have non-unique names but are in different subfolders. This overrides the `input` option. The format of the file is below. You should give the full path of each file. Sample names have to be unique and you should provide an entry for all the samples you want to run.

```
sample,filename1,filename2
S1,/folder/S1/cleaned_1.fastq.gz,/folder/S1/cleaned_2.fastq.gz
S2,/folder/S2/cleaned_1.fastq.gz,/folder/S2/cleaned_2.fastq.gz
S3,/folder/S3/cleaned_1.fastq.gz,/folder/S3/cleaned_2.fastq.gz
```

## Outputs

These files will be saved to the output folder when the workflow is finished.

```
raw.bcf - unfiltered output from bcftools mpileup, not overwritten by default
calls.vcf - unfiltered variant calls
filtered.vcf.gz - filtered vcf from all variant calls
snps.vcf.gz - snps only calls, used to make the core alignment
indels.vcf.gz - indels only, made from filtered calls
core.fa - fasta alignment from core snps, can be used to make a phylogeny
core.txt - text table of core snps
csq.tsv - consequence calls (if genbank provided)
csq_indels.tsv - consequence calls for indels
csq.matrix - matrix of consequence calls
snpdist.csv - comma separated distance matrix using snps
samples.csv - summary table of samples
RAxML_bipartitions.variants - ML tree if RAxML was used, optional
tree.newick - tree with SNPs branch lengths, if RAxMl used
```

## Use from Python

You can run a workflow from within Python by importing the snipgenie package and invoking the `WorkFlow` class. You need to provide the options in a dictionary with the same keywords as the command line. Notice in this example we are loading files from two folders.

```python
import snipgenie
args = {'threads':8, 'outdir': 'results', 'labelsep':'-',
        'input':['/my/folder/',
                 '/my/other/folder'],
        'reference': 'sequence.fa', 'overwrite':False}
W = snipgenie.app.WorkFlow(**args)
W.setup()
W.run()
```

## GUI

The package includes a desktop application with additional features like a fastq quality analysis, the ability to view alignments and tree viewing. It requires the installation of either PyQt5 or PySide2 if using the pip install. A windows installer for this application will be available separately.

<img src=doc/source/scr1.png width=600px>

You can view a short video on using the GUI [here](https://www.youtube.com/watch?v=ieljx3NiTqE).

## FAQ

_The run was stopped during execution, can it be resumed?_

Yes, by default the program won't overwrite intermediate files when re-run. So just run it again. Make sure there are no old tmp.****.bam files in the mapped folder if an alignment got interrupted. Variant calling can't be resumed though, so it will have to start again if interrupted.

_My sample files are not being parsed properly._

This may be because your sample names are unusual. The program extracts the unique sample names from the files by using the '_' symbol as delimeter. If your names differ you can supply a different delimeter with the labelsep option.

_I added new files and tried to re-run but it failed._

This is because the samples don't match the previous variant call output. You might see a `different number of samples` warning. By default the results of mpileup are not overwritten as this is the slowest step. You should first delete the file `raw.bcf` in the output folder and run again.

_I have more than 1000 samples and the bcftools mpileup step fails._

This is likely due to the limit on the number of files that can be opened at the same time. You can increase this limit on Linux using `ulimit -n 2000` or whatever value you need up to 9999. Note that for many samples this step could take several days to run.

## BTBGENIE

This is software is developed as part of the **BTBGENIE** project, a DAFM funded project for development of genomic epidemiology systems for tracking and eradicating Mycobacterium bovis in Ireland.
