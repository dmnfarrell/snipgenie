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

#location for sequences to align against
index_path = os.path.join(app.config_path, 'contam')
if not os.path.exists(index_path):
    os.makedirs(index_path, exist_ok=True)

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
        """Customise this and/or create_widgets"""

        if parent==None:
            return
        self.parent = parent
        self.outpath = os.path.join(self.parent.outputdir, 'contam')
        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath, exist_ok=True)

        self.refs = pd.DataFrame(columns=['name','description','filename','species'])
        self.ref_table = os.path.join(index_path,'contam_refs.csv')
        self.numreads = 100000
        self.create_widgets()
        self.fetch_defaults()
        if not os.path.exists(self.ref_table):
            self.find_files()
        self.load_table()
        #get previous results if present
        self.get_results()
        return

    def create_widgets(self):
        """Create widgets if GUI plugin"""

        self.main = QWidget()
        layout = self.layout = QVBoxLayout()
        layout.setAlignment(QtCore.Qt.AlignTop)
        self.main.setLayout(layout)

        self.tree = QTreeWidget()
        self.tree.setHeaderItem(QTreeWidgetItem(["name","species","desc","filename"]))
        self.tree.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tree.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self.show_tree_menu)
        layout.addWidget(self.tree)

        #t = self.result_table = tables.DataFrameTable(self.main)
        t = self.table_widget = tables.DataFrameWidget(self.main, toolbar=True)
        self.result_table = self.table_widget.table
        layout.addWidget(t)
        bw = self.create_buttons(self.main)
        layout.addWidget(bw)
        return

    def show_tree_menu(self, pos):
        """Show right cick tree menu"""

        item = self.tree.itemAt( pos )
        menu = QMenu(self.tree)
        propsAction = menu.addAction("Properties")
        deleteAction = menu.addAction("Delete")
        action = menu.exec_(self.tree.mapToGlobal(pos))
        if action == propsAction:
            self.edit_sequence_properties()
        elif action == deleteAction:
            self.remove_sequence()
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
        w.setText(str(self.numreads))
        vbox.addWidget(w)

        button = QPushButton("Set Folder")
        button.clicked.connect(self.set_folder)
        vbox.addWidget(button)
        button = QPushButton("Clean Files")
        button.clicked.connect(self.clean_files)
        vbox.addWidget(button)
        #button = QPushButton("Find NCBI Sequence")
        #button.clicked.connect(self.find_sequence_ncbi)
        #vbox.addWidget(button)
        button = QPushButton("Add sequence")
        button.clicked.connect(self.add_sequence)
        vbox.addWidget(button)
        button = QPushButton("Fetch NCBI Sequence")
        button.clicked.connect(self.fetch_sequence)
        vbox.addWidget(button)
        button = QPushButton("Run")
        button.setStyleSheet("background-color: red")
        button.clicked.connect(self.run)
        vbox.addWidget(button)
        return bw

    def update_folders(self):
        """Update output folders"""

        return

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

    def load_table(self):
        """load genomes table"""

        self.refs = pd.read_csv(self.ref_table, index_col=0)
        self.update_treelist()
        return

    def find_files(self):
        """Find sequence files in default folder and store them in table if missing"""

        for f in glob.glob(os.path.join(index_path,'*.fa')):
            #print (f)
            rec = SeqIO.read(f, "fasta")
            idx = os.path.basename(f)
            #print (rec)
            first, _, desc = rec.description.partition(" ")
            d = rec.description
            species = '_'.join(d.split()[1:3])
            self.refs.loc[rec.id] = [rec.id, rec.description, f, species]

        self.update_treelist()
        self.refs.to_csv(self.ref_table)
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
        """Get some default genome seqs from Genbank"""

        names = ['NC_002516','NZ_AP022609','CP089304','NZ_CP01680']
        for name in names:
            self.efetch(name)
        return

    def read_fasta_header(self, filename):
        """Fasta header only"""

        with open(filename, 'r') as f:
            x = f.readline().strip()[1:]
        return x.split(' ')

    def add_sequence(self):
        """Add a sequence from local fasta- may be more than one seq in fasta file"""

        filename, _ = QFileDialog.getOpenFileName(self.main, 'Add sequence', './',
                                        filter="Fasta Files(*.fa *.fasta);;All Files(*.*)")
        if not filename:
            return
        #header = self.read_fasta_header(filename)
        #print (header)
        recs = list(SeqIO.parse(filename, format='fasta'))
        rec = recs[0]

        d = rec.description
        species = '_'.join(d.split()[1:3])
        self.refs.loc[rec.id] = [rec.id, rec.description, filename, species]
        self.update_treelist()
        #print (self.refs)
        self.refs.to_csv(self.ref_table)
        return

    def remove_sequence(self):
        """Remove sequence"""

        item = self.tree.selectedItems()[0]
        row = self.tree.selectedIndexes()[0].row()
        name = item.text(0)
        self.refs = self.refs.drop(name)
        self.refs.to_csv(self.ref_table)
        self.tree.takeTopLevelItem(row)
        return

    def edit_sequence_properties(self):
        """Edit ref sequence"""

        item = self.tree.selectedItems()[0]
        name = item.text(0)
        g = self.refs.loc[name]

        opts = {'name':{'type':'entry','default':g.name},
                'desc':{'type':'entry','default':g.description},
                'species':{'type':'entry','default':g.species},
                }
        dlg = widgets.MultipleInputDialog(self.main, opts, title='Edit Genome Info',
                            width=500,height=150)
        dlg.exec_()
        if not dlg.accepted:
            return
        kwds = dlg.values
        self.refs = self.refs.drop(name)
        new = kwds['name']
        self.refs.loc[new] = [new,kwds['desc'],g.filename,kwds['species']]
        self.refs.to_csv(self.ref_table)
        #print (self.refs)
        self.update_treelist()
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

        def func(progress_callback):
            filename = os.path.join(index_path, accession+'.fa')
            if os.path.isfile(filename):
                print ('%s file present' %accession)
                return
            Entrez.email = "A.N.Other@example.com"
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            s = handle.read()
            with open(filename, "w") as out:
                out.write(s)
            handle.close()
            print ('downloaded to %s' %filename)
            rec = SeqIO.read(filename, "fasta")
            self.refs.loc[accession] = [rec.id, rec.description, filename, None]
        self.parent.run_threaded_process(func, self.fetch_completed)
        return

    def fetch_completed(self):

        self.update_treelist()
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

    def update_treelist(self):
        """Refresh entries in tree"""

        self.tree.clear()
        for i,r in self.refs.iterrows():
            #print (i,r)
            item = QTreeWidgetItem(self.tree)
            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(0, QtCore.Qt.Checked)
            item.setText(0, r.name)
            item.setText(1, r.species)
            item.setText(2, r.description)
            item.setText(3, r.filename)
        return

    def build_indexes(self):
        """Build index(es)"""

        def func(progress_callback):
            print ('building index..')
            #for a in genomes:
            for i,r in self.refs.iterrows():
                aligners.build_bwa_index(r.filename)
        self.parent.run_threaded_process(func, self.parent.processing_completed)
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

    def save_unmapped_reads(self):
        """Save unmapped reads to new files"""

        #cmd =
        return

    def run(self):
        """Run against selected genomes"""

        self.results = []
        def func(progress_callback):
            table = self.parent.fastq_table
            df = table.model.df.copy()
            rows = table.getSelectedRows()
            new = df.iloc[rows] #subset of sample table to run
            self.parent.opts.applyOptions()
            kwds = self.parent.opts.kwds
            threads = kwds['threads']
            self.numreads = int(self.readsentry.text())
            checked = self.get_checked()
            outdir = self.outpath
            #align to each index in turn
            for acc,g in self.refs.iterrows():
                if acc not in checked:
                    continue
                print (acc)
                aligners.build_bwa_index(g.filename, show_cmd=False, overwrite=False)
                for i,r in new.iterrows():
                    sample =  r['sample']
                    label = sample + '~~' + acc #label for bam
                    file1 = r.filename1
                    if 'filename2' in df.columns:
                        file2 = r.filename2
                    else:
                        file2 = None

                    #align subset reads
                    tmp1 = self.subset_reads(file1, outdir, self.numreads)
                    tmp2 = self.subset_reads(file2, outdir, self.numreads)
                    out = os.path.join(outdir,label+'.bam')
                    aligners.bwa_align(tmp1, tmp2, idx=g.filename, out=out,
                                        threads=threads, overwrite=False)
                    print (out)

        self.parent.run_threaded_process(func, self.run_completed)
        return

    def run_completed(self):
        self.parent.processing_completed()
        self.get_results()
        return

    def get_results(self):
        """Reuslts from mapping"""

        print ('getting results..')
        bamfiles = glob.glob(os.path.join(self.outpath,'*.bam'))
        #print (bamfiles)
        results = []
        checked = self.get_checked()
        for f in bamfiles:
            sample, idx = os.path.splitext(os.path.basename(f))[0].split('~~')
            if idx not in self.refs.index:
                continue
            r = self.refs.loc[idx]
            d = tools.samtools_flagstat(f)
            #total = tools.get_fastq_size(r.filename1)
            mapped = d['primary']
            perc = round(mapped/self.numreads*2*100,2) #assumes numreads unchanged since alignment!
            results.append([sample,idx,r.species,mapped,perc])
            #print (sample, idx)
            #print (d)

        df = self.results = pd.DataFrame(results, columns=['sample','ref','species','total','perc'])
        df = df[df.ref.isin(checked)]
        p = pd.pivot_table(df, index='sample',columns='ref',values='total')
        self.result_table.setDataFrame(p)
        return

    def save_data(self):
        """Return save data"""

        data = {}
        data['numreads'] = self.numreads
        data['checked'] = self.get_checked()
        return data

    def load_data(self, data):
        """Load any saved data from project. Run in constructor."""

        print (data)
        self.numreads = data['numreads']
        self.readsentry.setText(str(self.numreads))
        return

    def project_closed(self):
        """Run when parent project is closed"""

        
        return
