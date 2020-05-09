Using snpgenie
==============

The program includes both a command line and graphical interface. Both will produce the same results.

Command Line
------------

This will run the entire process based on a set of options given at the terminal::

  -h, --help            show this help message and exit
  -i FILE, --input FILE
                        input folder(s)
  -l FILE, --labels FILE
                        sample labels file, optional
  -r FILE, --reference FILE
                        reference genome filename
  -w, --overwrite       overwrite intermediate files
  -m, --trim            trim fastq files
  -q QUALITY, --quality QUALITY
                        trim quality
  -t THREADS, --threads THREADS
                        cpu threads to use
  -o FILE, --outdir FILE
                        Results folder
  -v, --version         Get version
  -s, --test            Do test run

Example::

  snpgenie -r reference.fa -g reference.gff -i data_files -t 8 -o results

From Python
-----------

You can run a workflow from within Python::

  from sngenie import app
  args = {'threads':8, 'outdir': 'results', 'labelsep':'-',
          'input':['/my/folder/',
                   '/my/other/folder'],
          'reference': None, 'overwrite':False}
  W = app.WorkFlow(**args)
  st = W.setup()
  W.run()

Desktop Application
-------------------

This interactive tool is designed for those not comfortable with the command line and includes some additional features such as visualization of fastq qualities.
