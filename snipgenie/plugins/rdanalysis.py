"""
    RD plugin for snipgenie
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
import pickle, gzip, subprocess
import random
from collections import OrderedDict
from snipgenie.qt import *
import pandas as pd
from Bio import Phylo, SeqIO
from Bio import Entrez
from snipgenie import app, widgets, tables, aligners, tools, rdiff
from snipgenie.plugin import Plugin

class RegionDiffPlugin(Plugin):
    """Region Diff plugin for SNiPgenie"""

    #uncomment capabilities list to appear in menu
    capabilities = ['gui','docked']
    requires = ['']
    menuentry = 'MTBC RD analysis'
    name = 'MTBC RD'
    iconfile = 'rd.svg'
    side = 'right' #dock location to load plugin

    def __init__(self, parent=None):
        """Customise this and/or create_widgets"""

        if parent==None:
            return
        self.parent = parent
        self.outpath = os.path.join(self.parent.outputdir, 'rd_aligned')
        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath, exist_ok=True)

        self.numreads = 100000
        self.create_widgets()
        #self.load_table()
        #get previous results if present
        #self.get_results()
        return

    def create_widgets(self):
        """Create widgets if GUI plugin"""

        self.main = QWidget()
        layout = self.layout = QVBoxLayout()
        layout.setAlignment(QtCore.Qt.AlignTop)
        self.main.setLayout(layout)

        t = self.table_widget = tables.DataFrameWidget(self.main, toolbar=True)
        self.result_table = self.table_widget.table
        layout.addWidget(t)
        bw = self.create_buttons(self.main)
        layout.addWidget(bw)
        return

    def create_buttons(self, parent):

        bw = QWidget(parent)
        vbox = QVBoxLayout(bw)

        w = QLabel('Read limit:')
        vbox.addWidget(w)
        self.readsentry = w = QLineEdit()
        w.setText(str(self.numreads))
        vbox.addWidget(w)
        button = QPushButton("Set Folder")
        button.clicked.connect(self.set_folder)
        vbox.addWidget(button)
        button = QPushButton("Clean Files")
        button.clicked.connect(self.clean_files)
        vbox.addWidget(button)
        button = QPushButton("Run")
        button.clicked.connect(self.run)
        vbox.addWidget(button)
        return bw

    def clean_files(self):
        """Clean aligned files"""

        reply = QMessageBox.question(self.main, 'Confirm',
                                "This will clean all aligned files.\nAre you sure?",
                                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.No:
            return False
        files = glob.glob(os.path.join(self.outpath,'*.*'))
        for f in files:
            os.remove(f)
        return

    def set_folder(self):
        """Set test folder"""

        path = QFileDialog.getExistingDirectory(self.main, 'Set Folder', './')
        if not path:
            return
        self.outpath = path
        self.update_folders()
        return


    def subset_reads(self, filename, path='/tmp', numreads=10000, overwrite=False):
        """Subset of reads"""

        from Bio import SeqIO, bgzf
        from gzip import open as gzopen
        new = os.path.join(path,os.path.basename(filename))
        if os.path.exists(new) and overwrite == False:
            return new
        print ('subsetting %s reads' %numreads)
        recs = SeqIO.parse(gzopen(filename, "rt"), "fastq")
        i=0
        with bgzf.BgzfWriter(new, "wb") as outgz:
            for rec in recs:
                if i>numreads:
                    break
                SeqIO.write(sequences=rec, handle=outgz, format="fastq")
                i+=1
        return new

    def run(self):
        """Run against selected genomes"""

        self.results = []
        def func(progress_callback):
            table = self.parent.fastq_table
            #df = table.model.df.copy()
            #rows = table.getSelectedRows()
            #new = df.iloc[rows] #subset of sample table to run
            self.parent.opts.applyOptions()
            kwds = self.parent.opts.kwds
            threads = kwds['threads']

            rdiff.create_rd_index()
            data = self.parent.get_selected()
            if data is None or len(data) == 0:
                return
            out = self.outpath
            res = rdiff.run_samples(data, out, threads=kwds['threads'])

            X = rdiff.get_matrix(res, cutoff=0.15)
            X['species'] = X.apply(rdiff.apply_rules,1)
            self.rd_result = X

        self.parent.run_threaded_process(func, self.run_completed)
        return

    def run_completed(self):
        self.parent.processing_completed()
        self.get_results()
        return

    def get_results(self):
        """Reuslts from mapping"""

        print ('getting results..')

        self.result_table.setDataFrame(self.rd_result)
        return

    def save_data(self):
        """Return save data"""

        data = {}
        return data

    def load_data(self, data):
        """Load any saved data from project.
        Run when plugin is initially launched."""

        return
