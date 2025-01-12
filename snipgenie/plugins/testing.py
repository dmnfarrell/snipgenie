"""
    Testing plugin for snipgenie
    Created October 2022
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

import sys,os,platform,time,tempfile,glob
import pickle, gzip
import random
from collections import OrderedDict
from snipgenie.qt import *
import pandas as pd
from Bio import Phylo
from snipgenie import app, widgets, tables
from snipgenie.plugin import Plugin
from snipgenie import simulate

newick = os.path.join(app.module_path, 'testing', 'sim.newick')

class TestingPlugin(Plugin):
    """Testing plugin for SNiPgenie"""

    #uncomment capabilities list to appear in menu
    capabilities = ['gui','docked']
    requires = ['']
    menuentry = 'Testing'
    name = 'Testing'
    iconfile = 'tests.svg'
    side = 'right' #dock location to load plugin
    enabled = False
    
    def __init__(self, parent=None, table=None):
        """Customise this and/or doFrame for your widgets"""

        if parent==None:
            return
        self.parent = parent
        self.table = table
        self.outpath =  tempfile.gettempdir()
        self.refs = ['Mbovis-AF212297','Sars-Cov-2']
        self.update_folders()
        self.create_widgets()
        return

    def create_widgets(self):
        """Create widgets if GUI plugin"""

        self.main = QWidget()
        layout = self.layout = QVBoxLayout()
        layout.setAlignment(QtCore.Qt.AlignTop)
        self.main.setLayout(layout)
        #add buttons
        bw = self.create_buttons(self.main)
        layout.addWidget(bw)
        return

    def create_buttons(self, parent):

        bw = QWidget(parent)
        vbox = QVBoxLayout(bw)
        #button = QPushButton("Download Test Data")
        #button.clicked.connect(self.download)
        #vbox.addWidget(button)
        w = QLabel('Reference:')
        vbox.addWidget(w)
        self.refbox = w = QComboBox()
        w.addItems(self.refs)
        w.setCurrentIndex(0)
        vbox.addWidget(w)
        w = QLabel('Reads:')
        vbox.addWidget(w)
        self.readsentry = w = QLineEdit()
        w.setText(str(100000))
        vbox.addWidget(w)

        button = QPushButton("Set Folder")
        button.clicked.connect(self.set_folder)
        vbox.addWidget(button)
        button = QPushButton("Clear Files")
        button.clicked.connect(self.clear_files)
        vbox.addWidget(button)
        button = QPushButton("Simulate Reads")
        button.clicked.connect(self.simulate_reads)
        vbox.addWidget(button)
        button = QPushButton("Load test files")
        button.clicked.connect(self.load_files)
        vbox.addWidget(button)
        button = QPushButton("Compare Trees")
        button.clicked.connect(self.compare_tree)
        vbox.addWidget(button)
        return bw

    def update_folders(self):
        """Update output folders"""

        self.testfasta = os.path.join(self.outpath,'phastsim_output1.fasta')
        self.testdir = os.path.join(self.outpath,'sim_fastq')
        self.resultsdir = os.path.join(self.outpath,'sim_results')
        return

    def download(self):
        """Fetch test data"""

        url = 'https://zenodo.org/record/5179838/files/Artificial_reads.zip'
        return

    def set_folder(self):
        """Set test folder"""

        path = QFileDialog.getExistingDirectory(self.main, 'Set Folder', './')
        if not path:
            return
        self.outpath = path
        self.update_folders()
        return

    def clear_files(self):
        """Delete test files"""

        files = glob.glob(os.path.join(self.testdir,'*'))
        files += glob.glob(os.path.join(self.resultsdir,'*'))
        files.extend(glob.glob(os.path.join(self.outpath,'phastsim_output*')))
        for f in files:
            os.remove(f)
        return

    def simulate_reads(self):
        """Simulate artificial reads"""

        name = self.refbox.currentText()
        #add ref genome
        ref = app.preset_genomes[name]['sequence']
        print ('using %s' %ref)
        reads = int(self.readsentry.text())
        def func(progress_callback):
            print ('simulate genomes..')
            if not os.path.exists(self.testfasta):
                simulate.run_phastsim(self.outpath, ref, newick)
            print ('making fastq files..')
            outpath = self.testdir
            os.makedirs(outpath, exist_ok=True)
            infile = self.testfasta
            simulate.generate_fastqs(infile, outpath, reads=reads, overwrite=False)

        self.parent.run_threaded_process(func, self.parent.processing_completed)
        return

    def load_files(self):
        """Load the test files into a new project"""

        parent = self.parent
        os.makedirs(self.resultsdir, exist_ok=True)
        parent.outputdir = self.resultsdir
        parent.labelsep = '.'
        parent.update_labels()
        filenames = app.get_files_from_paths(self.testdir)
        #print (filenames)
        parent.opts.setWidgetValue('labelsep', '.')
        parent.load_fastq_table(filenames)
        self.create_meta_data()

        #add meta data
        df=parent.fastq_table.model.df
        df = df.merge(self.meta,left_index=True,right_index=True)
        parent.fastq_table.setDataFrame(df)

        ref = self.refbox.currentText()
        #add ref genome
        gm = app.preset_genomes[ref]
        parent.load_preset_genome(gm['sequence'], gm['gb'], gm['mask'], ask=False)
        return

    def compare_tree(self):
        """Compare result from workflow with test tree"""

        tree1 = Phylo.read(newick, "newick")
        Phylo.draw_ascii(tree1)
        tree2 = Phylo.read(parent.treefile, "newick")
        term_names1 = [term.name for term in tree1.get_terminals()]
        term_names2 = [term.name for term in tree2.get_terminals()]
        # false if terminals are not the same
        if set(term_names1) != set(term_names2):
            return False
        # true if _BitStrings are the same
        if _bitstrs(tree1) == _bitstrs(tree2):
            return True
        else:
            return False

        return

    def quit(self, evt=None):
        """Override this to handle pane closing"""

        self.main.close()
        return

    def about(self):
        """About this plugin"""

        txt = "This plugin implements ...\n"+\
               "version: %s" %self.version
        return txt
