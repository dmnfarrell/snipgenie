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
import pickle, gzip, subprocess
import random
from collections import OrderedDict
from snipgenie.qt import *
import pandas as pd
from Bio import Phylo, SeqIO
from Bio import Entrez
from snipgenie import app, widgets, tables, aligners, tools
from snipgenie.plugin import Plugin

index_path = os.path.join(app.config_path, 'contam')

class Genome(object):
    def __init__(self, name, description=None, filename=None):

        self.name = name
        self.filename = filename
        self.description = description
        return

class ContaminationCheckerPlugin(Plugin):
    """Contam checker plugin for SNiPgenie"""

    #uncomment capabilities list to appear in menu
    capabilities = ['gui']
    requires = ['']
    menuentry = 'Contamination Check'
    name = 'Contamination Check'
    side = 'right' #dock location to load plugin

    def __init__(self, parent=None, table=None):
        """Customise this and/or doFrame for your widgets"""

        if parent==None:
            return
        self.parent = parent
        self.table = table
        self.outpath = os.path.join(tempfile.gettempdir(), 'contam')
        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath, exist_ok=True)

        self.genomes = {}
        self.create_widgets()
        self.fetch_defaults()
        self.find_files()
        return

    def create_widgets(self):
        """Create widgets if GUI plugin"""

        self.main = QWidget()
        layout = self.layout = QVBoxLayout()
        layout.setAlignment(QtCore.Qt.AlignTop)
        self.main.setLayout(layout)
        #add buttons
        self.tree = QTreeWidget()
        self.tree.setHeaderItem(QTreeWidgetItem(["name","desc"]))
        layout.addWidget(self.tree)
        bw = self.create_buttons(self.main)
        layout.addWidget(bw)
        #self.createMenu(self.main)
        return

    def createMenu(self, parent):
        """Main menu"""

        self.menubar = QMenuBar(parent)
        self.entrez_menu = QMenu('Entrez', parent)
        self.menubar.addMenu(self.entrez_menu)
        self.entrez_menu.addAction('Default', self.find_sequence_ncbi)

        return

    def create_buttons(self, parent):

        bw = QWidget(parent)
        vbox = QVBoxLayout(bw)

        w = QLabel('Read limit:')
        vbox.addWidget(w)
        self.readsentry = w = QLineEdit()
        w.setText(str(100000))
        vbox.addWidget(w)

        button = QPushButton("Set Folder")
        button.clicked.connect(self.set_folder)
        vbox.addWidget(button)
        button = QPushButton("Clean Files")
        button.clicked.connect(self.clean_files)
        vbox.addWidget(button)
        button = QPushButton("Find Genome")
        button.clicked.connect(self.find_sequence_ncbi)
        vbox.addWidget(button)
        button = QPushButton("Fetch Sequence")
        button.clicked.connect(self.fetch_sequence)
        vbox.addWidget(button)
        button = QPushButton("Run")
        button.clicked.connect(self.run)
        vbox.addWidget(button)
        return bw

    def clean_files(self):
        """Clean"""

        return

    def find_files(self):
        """Find genome files"""

        for f in glob.glob(os.path.join(index_path,'*.fa')):
            print (f)
            rec = SeqIO.read(f, "fasta")
            idx = os.path.basename(f)
            #print (rec)
            self.genomes[rec.id] = Genome(rec.id, rec.description, f)
        print (self.genomes)
        self.update_sequences()
        return

    def update_folders(self):
        """Update output folders"""

        self.mapped = os.path.join(self.outpath, 'mapped')
        return

    def set_folder(self):
        """Set test folder"""

        path = QFileDialog.getExistingDirectory(self.main, 'Set Folder', './')
        if not path:
            return
        self.outpath = path
        self.update_folders()
        return

    def find_sequence_ncbi(self):
        """Search for sequence using esearch"""

        opts = {'keywords':{'type':'entry','default':''},
                }
        dlg = widgets.MultipleInputDialog(self.main, opts, title='Find Sequence',
                            width=250,height=150)
        dlg.exec_()
        if not dlg.accepted:
            return
        kwds = dlg.values
        acc = kwds['accession']
        keyword = kwds['keyword']
        #self.fetch_sequence(keyword)
        self.esearch(keyword)
        return

    def esearch(self, keyword=""):
        """entrez search"""

        def func(progress_callback):
            Entrez.email = "A.N.Other@example.com"
            handle = Entrez.esearch(db="nucleotide", term=keyword)
            record = Entrez.read(handle)
            print (record["IdList"])
        self.parent.run_threaded_process(func, self.parent.processing_completed)
        return

    def search_completed(self):

        return

    def fetch_defaults(self):
        """get default genome seqs"""

        names = ['NC_002516','NZ_AP022609','NZ_SILH01000002','CP089304','NZ_CP01680']
        for name in names:
            self.efetch(name)
        return

    def fetch_sequence(self):

        opts = {'accession':{'type':'entry','default':''},
                }
        dlg = widgets.MultipleInputDialog(self.main, opts, title='Find Sequence',
                            width=250,height=150)
        dlg.exec_()
        if not dlg.accepted:
            return
        kwds = dlg.values
        self.efetch(kwds['accession'])
        return

    def efetch(self, accession):
        """Fetch a fasta sequence from genbank nucl DB"""

        print (accession)
        def func(progress_callback):
            filename = os.path.join(index_path, accession+'.fa')
            if os.path.isfile(filename):
                return
            Entrez.email = "A.N.Other@example.com"
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            s = handle.read()
            with open(filename, "w") as out:
                out.write(s)
            handle.close()
            print ('downloaded to %s' %filename)
            rec = SeqIO.read(filename, "fasta")
            #print (rec.description)
            self.genomes[accession] = Genome(rec.id, rec.description, filename)
        self.parent.run_threaded_process(func, self.fetch_completed)
        return

    def fetch_completed(self):

        self.update_sequences()
        self.parent.processing_completed()
        return

    def get_checked(self):

        names=[]
        for i in range(self.tree.topLevelItemCount()):
            item = self.tree.topLevelItem(i)
            if item.checkState(0) == QtCore.Qt.CheckState.Checked:
                names.append(item.text(0))
        return names

    def get_items(self):

        names=[]
        for i in range(self.tree.topLevelItemCount()):
            item = self.tree.topLevelItem(i)
        return names

    def update_sequences(self):
        """Add entry to tree"""

        names = self.get_items()
        for i in self.genomes:
            if i in names:
                continue
            g = self.genomes[i]
            item = QTreeWidgetItem(self.tree)
            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(0, QtCore.Qt.Checked)
            item.setText(0, i)
            item.setText(1, g.description)

        return

    def build_indexes(self):
        """Build index(es)"""

        def func(progress_callback):
            print ('building index..')
            for a in genomes:
                g=genomes[a]
                aligners.build_bwa_index(g.filename)
        self.parent.run_threaded_process(func, self.parent.processing_completed)
        return

    def run(self):
        """Run against genomes"""

        table = self.parent.fastq_table
        df = table.model.df.copy()
        rows = table.getSelectedRows()
        new = df.iloc[rows]
        self.mapped = os.path.join(self.outpath, 'mapped')
        threads = 8
        def func(progress_callback):
            for idx in self.genomes:
                print (idx)
                g = self.genomes[idx]
                aligners.build_bwa_index(g.filename)
                self.df = app.align_reads(new, idx=g.filename, outdir=self.mapped, overwrite=False)
                print(self.df)
        self.parent.run_threaded_process(func, self.parent.processing_completed)
        #print (self.df)
        #app.mapping_stats(samples)
        return

    def about(self):
        """About this plugin"""

        txt = "This plugin is a contamination checker...\n"+\
               "version: %s" %self.version
        return txt
