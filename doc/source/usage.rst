Usage
=====

The program includes both a command line and graphical interface. Both will produce the same results.

Command Line
------------

This will run the entire process based on a set of options given at the terminal::

  -h, --help            show this help message and exit
  -i FILE, --input FILE
                        input folder(s)
  -e LABELSEP, --labelsep LABELSEP
                        symbol to split the sample labels on
  -r FILE, --reference FILE
                        reference genome filename
  -S SPECIES, --species SPECIES
                        set the species reference genome, overrides -r
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
                        sequencing platform, change to ont if using oxford nanopore
  -a ALIGNER, --aligner ALIGNER
                        aligner to use, bwa, subread, bowtie or minimap2
  -b, --buildtree       whether to build a phylogenetic tree, requires RaXML
  -N BOOTSTRAPS, --bootstraps BOOTSTRAPS
                        number of bootstraps to build tree
  -o FILE, --outdir FILE
                        Results folder
  -q, --qc              Get version
  -d, --dummy           Check samples but don't run
  -x, --test            Test run
  -v, --version         Get version

Examples::

Call with your own reference fasta file::

  snipgenie -r reference.fa -i data_files -o results

Use an in built species genome as reference. This will also supply an annotation file. The current options are Mbovis-AF212297, MTB-H37Rv, MAP-K10, M.smegmatis-MC2155::

  snipgenie -S Mbovis-AF212297 -i data_files -o results

Provide more than one folder::

  snipgenie -r reference.fa -i data_files1 -i data_files2 -o results

Provide an annotation (genbank format) for consequence calling::

  snipgenie -r reference.fa -g reference.gb -i data_files -o results

Add your own filters and provide threads::

  snipgenie -r reference.fa -i data_files -t 8 -o results` \
   -f 'QUAL>=40 && INFO/DP>=20 && MQ>40'

Aligners
++++++++

You can use any one of the following aligners: bwa, subread, bowtie or minimap2. These should be present on your system, unless using the Windows version. Note that for oxford nanopore reads you should use minimap2 and specify the platform as 'ont'.

Mask file
+++++++++

You can selectively mask snp sites such as those contained in transposons or repetitive regions from being included in the output. You need to provide a bed file with the following columns: chromosome name, start and end coordinates of the regions. There is currently a built-in mask file used for M.bovis and of you select this genome as reference using the --species option it will be used automatically.
Example::

  LT708304.1 	 105359 	 106751
  LT708304.1 	 131419 	 132910
  LT708304.1 	 149570 	 151187
  LT708304.1 	 306201 	 307872

Inputs
++++++

You can provide single folder with all the files in one place or multiple folders. Folders are searched recursively for inputs with extensions *.f*q.gz. So be careful you don't have files in the folders you don't want included. The following file structure will load both sets of files if you provide the parent folder as input. You can also provide multiple separate folders using -i as shown above.

For example if you provide -i data with the following structure::

  data/
  ├── ERR1588781
  │   ├── ERR1588781_1.fq.gz
  │   └── ERR1588781_2.fq.gz
  └── ERR1588785
      ├── ERR1588785_1.fastq.gz
      └── ERR1588785_2.fastq.gz

Filenames are parsed and a sample name is extracted for each pair (if paired end). This is simply done by splitting on the _ symbol. So a file called /path/13-11594_S85_L001-4_R1_001.fastq.gz will be given a sample name 13-11594. As long as the sample names are unique this is ok. If you had a file names like A_2_L001-4_R1_001, A_3_L001-4_R1_001 you should split on '-' instead. You can specify this in the labelsep option. The workflow won't run unless sample names are unique.

Outputs
+++++++

These files will be saved to the output folder when the workflow is finished::

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

Use from Python
---------------

Run the workflow
++++++++++++++++

You can run a workflow from within Python by importing the snipgenie package and invoking the WorkFlow class. You need to provide the options in a dictionary with the same keywords as the command line. Notice in this example we are loading files from two folders.
::

  from snipgenie import app
  args = {'threads':8, 'outdir': 'results', 'labelsep':'-',
          'input':['/my/folder/',
                   '/my/other/folder'],
          'reference': None, 'overwrite':False}
  W = app.WorkFlow(**args)
  st = W.setup()
  W.run()

Run previously aligned files
++++++++++++++++++++++++++++

This will run the remainder of the variant calling from sets of bam files. This is useful if you want to put together a set of precious runs without aligning them all again::

  from snipgenie import app
  #get a list of bam files
  bams = glob.glob('/my/folder/*.bam')
  ref = app.mbovis_genome # use your own reference here
  app.run_bamfiles(bams, ref)
