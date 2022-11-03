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

    def get_species(self):
        """Get species from desc"""

        d = self.description
        self.species = '_'.join(d.split()[:2])
        return

    def __repr__(self):
        return 'genome object %s' %self.__dict__

class ContaminationCheckerPlugin(Plugin):
    """Contam checker plugin for SNiPgenie"""

    #uncomment capabilities list to appear in menu
    capabilities = ['gui']
    requires = ['']
    menuentry = 'Contamination Check'
    name = 'Contamination Check'
    iconfile = 'contam.png'
    side = 'right' #dock location to load plugin

    def __init__(self, parent=None):
        """Customise this and/or doFrame for your widgets"""

        if parent==None:
            return
        self.parent = parent
        self.outpath = os.path.join(self.parent.outputdir, 'contam')
        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath, exist_ok=True)
        self.mapped = os.path.join(self.outpath, 'mapped')
        if not os.path.exists(self.mapped):
            os.makedirs(self.mapped, exist_ok=True)

        self.genomes = {}
        self.create_widgets()
        self.fetch_defaults()
        self.find_files()
        #get previous results if present
        self.get_results()
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
        t = self.result_table = tables.DataFrameTable(self.main)
        #t.resize(20,5)
        layout.addWidget(t)
        bw = self.create_buttons(self.main)
        layout.addWidget(bw)
        #self.createMenu(self.main)
        return

    def create_menu(self, parent):
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
        #button = QPushButton("Find Genome")
        #button.clicked.connect(self.find_sequence_ncbi)
        #vbox.addWidget(button)
        button = QPushButton("Fetch Sequence")
        button.clicked.connect(self.fetch_sequence)
        vbox.addWidget(button)
        button = QPushButton("Run")
        button.clicked.connect(self.run)
        vbox.addWidget(button)
        return bw

    def update_folders(self):
        """Update output folders"""

        self.mapped = os.path.join(self.outpath, 'mapped')
        return

    def clean_files(self):
        """Clean aligned files"""

        files = glob.glob(os.path.join(self.mapped,'*.bam'))
        for f in files:
            os.remove(f)
        return

    def find_files(self):
        """Find genome files and store them"""

        for f in glob.glob(os.path.join(index_path,'*.fa')):
            #print (f)
            rec = SeqIO.read(f, "fasta")
            idx = os.path.basename(f)
            #print (rec)
            first, _, desc = rec.description.partition(" ")
            g = self.genomes[rec.id] = Genome(rec.id, desc, f)
            g.get_species()
        #print (self.genomes)
        self.update_sequences()
        return

    def set_folder(self):
        """Set test folder"""

        path = QFileDialog.getExistingDirectory(self.main, 'Set Folder', './')
        if not path:
            return
        self.outpath = path
        self.update_folders()
        return

    def fetch_defaults(self):
        """get default genome seqs"""

        names = ['NC_002516','NZ_AP022609','CP089304','NZ_CP01680']
        for name in names:
            self.efetch(name)
        return

    def fetch_sequence(self):
        """Fetch a fasta sequence from genbank nucl DB"""

        opts = {'accession':{'type':'entry','default':''},
                }
        dlg = widgets.MultipleInputDialog(self.main, opts, title='Enter NCBI Accession',
                            width=250,height=150)
        dlg.exec_()
        if not dlg.accepted:
            return
        kwds = dlg.values
        self.efetch(kwds['accession'])
        return

    def efetch(self, accession):
        """Fetch a fasta sequence from entrez using efetch"""

        print (accession)
        def func(progress_callback):
            filename = os.path.join(index_path, accession+'.fa')
            if os.path.isfile(filename):
                print ('file present')
                return
            Entrez.email = "A.N.Other@example.com"
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            s = handle.read()
            with open(filename, "w") as out:
                out.write(s)
            handle.close()
            print ('downloaded to %s' %filename)
            rec = SeqIO.read(filename, "fasta")
            self.genomes[accession] = Genome(rec.id, rec.description, filename)
        self.parent.run_threaded_process(func, self.fetch_completed)
        return

    def fetch_completed(self):

        self.update_sequences()
        self.parent.processing_completed()
        return

    def find_sequence_ncbi(self):
        """Search for sequence using esearch"""

        opts = {'terms':{'type':'entry','default':'','label':'Search Terms'},
                }
        dlg = widgets.MultipleInputDialog(self.main, opts, title='Find Sequence',
                            width=250,height=150)
        dlg.exec_()
        if not dlg.accepted:
            return
        kwds = dlg.values
        terms = kwds['terms']
        self.esearch(terms)
        return

    def esearch(self, terms):
        """entrez search"""

        def func(progress_callback):
            Entrez.email = "A.N.Other@example.com"
            handle = Entrez.esearch(db="nucleotide", term=terms)
            record = Entrez.read(handle)
            print (record["IdList"])
        self.parent.run_threaded_process(func, self.search_completed)
        return

    def search_completed(self):
        self.parent.processing_completed()
        return

    def get_checked(self):
        """Get checked sequence names"""

        names=[]
        for i in range(self.tree.topLevelItemCount()):
            item = self.tree.topLevelItem(i)
            if item.checkState(0) == QtCore.Qt.CheckState.Checked:
                names.append(item.text(0))
        return names

    def get_items(self):

        names=[]
        #print (self.tree.topLevelItemCount())
        for i in range(self.tree.topLevelItemCount()):
            item = self.tree.topLevelItem(i)
            names.append(item.text(0))
        return names

    def update_sequences(self):
        """Refresh entries in tree"""

        for i in self.genomes:
            names = self.get_items()
            if i in names:
                continue
            #print (names, i)
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
        """Run against selected genomes"""

        self.results = []
        def func(progress_callback):
            table = self.parent.fastq_table
            df = table.model.df.copy()
            rows = table.getSelectedRows()
            new = df.iloc[rows]
            threads = 8
            names = self.get_checked()
            outdir = self.mapped
            #align to each index in turn
            for acc in self.genomes:
                print (acc)
                g = self.genomes[acc]
                aligners.build_bwa_index(g.filename, show_cmd=False)
                for i,r in new.iterrows():
                    sample =  r['sample']
                    label = sample + '~~' + acc #label for bam
                    file1 = r.filename1
                    if 'filename2' in df.columns:
                        file2 = r.filename2
                    else:
                        file2 = None

                    out = os.path.join(outdir,label+'.bam')
                    aligners.bwa_align(file1, file2, idx=g.filename, out=out, overwrite=False)
                    print (out)

        self.parent.run_threaded_process(func, self.run_completed)
        return

    def run_completed(self):
        self.parent.processing_completed()
        self.get_results()
        return

    def get_results(self):
        """Reuslts from mapping"""

        bamfiles = glob.glob(os.path.join(self.mapped,'*.bam'))
        #print (bamfiles)
        results = []
        for f in bamfiles:
            sample, idx = os.path.splitext(os.path.basename(f))[0].split('~~')
            d = tools.samtools_flagstat(f)
            total = d['total']
            g = self.genomes[idx]
            results.append([sample,idx,g.species,total])
            print (sample, idx)
            #print (d)

        df = self.results = pd.DataFrame(results, columns=['sample','ref','species','total'])
        print (df)
        p = pd.pivot_table(df, index='sample',columns='species',values='total')
        self.result_table.setDataFrame(p)
        return
