"""
    Species 16S genes check plugin for snipgenie
    Created April 2024
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

import sys,os,platform,time
import shutil,tempfile,glob
import pickle, gzip, subprocess
import random
from collections import OrderedDict
from snipgenie.qt import *
import pandas as pd
from Bio import Phylo, SeqIO

from snipgenie import app, widgets, tables, tools, trees
from snipgenie.plugin import Plugin

#location for sequences to align against
ref_path = os.path.join(app.config_path, 'species_16S')
if not os.path.exists(ref_path):
    os.makedirs(ref_path, exist_ok=True)
ref_file = os.path.join(ref_path,'16S_ncbi.fa')

def get_blast_coverage(bl, fasta):
    """Get alignment coverage of blast results from original sequence lengths"""

    df=tools.fasta_to_dataframe('16S_ncbi.fa')
    df=df.rename(columns={'length':'sslen'})
    bl=bl.merge(df,left_on='sseqid',right_on='name',how='left')
    bl['perc_cov'] = bl.apply(lambda x: round(x.length/x.sslen*100,2),1)
    return bl

def blast_16S(filename, hits=100, pident=99.5):
    """blast 16S genes"""

    bl=tools.blast_fasta('16S_ncbi.fa',filename,maxseqs=hits)
    bl = get_blast_coverage(bl, '16S_ncbi.fa')
    #bl['genus'] = bl.stitle.apply(lambda x: x.split()[1])
    bl['species'] = bl.stitle.apply(lambda x: ' '.join(x.split()[1:3]))
    cols = ['sseqid','sslen','length','perc_cov','pident','stitle','species']
    bl = bl[(bl.pident>=pident) & (bl.perc_cov>=80)].sort_values('pident',ascending=False)
    return bl

def extract_sequences_by_ids(input_fasta, output_fasta, ids_to_extract):
    """Extract sequences from fasta file with given ids"""

    sequences = SeqIO.parse(input_fasta, "fasta")
    # Filter sequences that match the given IDs
    filtered = (seq for seq in sequences if seq.id in ids_to_extract)
    SeqIO.write(filtered, output_fasta, "fasta")
    return

def append_sequences_to_fasta(fasta_file, seqs):
    """Append SeqRecords to a FASTA file, overwriting the old file."""

    existing_seqs = list(SeqIO.parse(fasta_file, "fasta"))
    if type(seqs) is not list:
        seqs = [seqs]
    existing_seqs.extend(seqs)
    #with open(fasta_file, "w") as output_handle:
    SeqIO.write(existing_seqs, fasta_file, "fasta")
    return

#extract_sequences_by_ids('16S_ncbi.fa', 'temp.fa',  ['NR_115676.1','NR_179517.1'])

def get_tree(fasta_file):
    """get phylo tree from fasta_file with fasttree"""

    out = 'temp.newick'
    cmd=f'mafft {fasta_file} > temp.aln'
    tmp = subprocess.check_output(cmd, shell=True, stderr=subprocess.PIPE)
    cmd=f'fasttree temp.aln > {out}'
    tmp = subprocess.check_output(cmd, shell=True, stderr=subprocess.PIPE)
    return out

class SpeciesCheckerPlugin(Plugin):
    """16S species checker plugin for SNiPgenie"""

    #uncomment capabilities list to appear in menu
    capabilities = ['gui','docked']
    requires = ['']
    menuentry = 'Species Check'
    name = '16S Species Check'
    iconfile = '16S.svg'
    side = 'right' #dock location to load plugin

    def __init__(self, parent=None):
        """Customise this and/or create_widgets"""

        if parent==None:
            return
        self.parent = parent
        if self.parent.outputdir == None:
            print ('no output dir set in project!')
            return
        self.outpath = os.path.join(self.parent.outputdir, 'species_16S')
        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath, exist_ok=True)

        self.numreads = 100000
        self.mismatches = 4
        self.create_widgets()
        self.fetch_sequence()
        return

    def fetch_sequence(self):
        """Get 16S sequences"""

        if os.path.exists(ref_file):
            return
        print ('fetching 16S sequence data..')
        url = 'https://github.com/dmnfarrell/snipgenie/raw/master/extra/16S_ncbi.fa.gz'
        infile = ref_file+'.gz'
        import urllib.request
        urllib.request.urlretrieve(url, infile)
        tools.gunzip(infile, ref_file)
        tools.make_blast_database(ref_file)
        return

    def create_widgets(self):
        """Create widgets if GUI plugin"""

        self.main = QWidget()
        layout = self.layout = QVBoxLayout()
        layout.setAlignment(QtCore.Qt.AlignTop)
        self.main.setLayout(layout)

        t = self.table_widget = tables.DataFrameWidget(self.main, app=self.parent, toolbar=True)
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
        self.readsentry = w = QSpinBox()
        w.setRange(10000,1e6)
        #w.setMaximum(max)
        #w.setMinimum(min)
        w.setSingleStep(10000)
        w.setValue(self.numreads)
        vbox.addWidget(w)
        w = QLabel('Use percentage:')
        vbox.addWidget(w)
        self.percbox = w = QCheckBox()
        w.setChecked(0)
        vbox.addWidget(w)
        #w = QLabel('Mismatches:')
        vbox.addWidget(w)
        #self.mismatchentry = w = QSpinBox()
        #w.setValue(self.mismatches)
        #vbox.addWidget(w)
        button = QPushButton("Clean Files")
        button.clicked.connect(self.clean_files)
        vbox.addWidget(button)
        button = QPushButton("Run")
        button.setStyleSheet("background-color: red")
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
            outdir = self.outpath
            for i,r in new.iterrows():
                if not os.path.exists(r.assembly):
                    continue
                bl = blast_16S(f{r.assembly},pident=90,hits=200)
                print (bl)

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
        results = []

        #df = self.results = pd.DataFrame(results, columns=['sample','ref','species','total','perc'])
  

        #p = pd.pivot_table(df, index='sample',columns='species',values=valcol)
        #self.result_table.setDataFrame(p)
        return

    def save_data(self):
        """Return save data"""

        data = {}    
        #data['checked'] = self.get_checked()
        return data

    def load_data(self, data):
        """Load any saved data from project.
        Run when plugin is initially launched."""

        self.numreads = data['numreads']
        self.readsentry.setValue(self.numreads)
        checked = data['checked']
        #print (checked)
        for i in range(self.tree.topLevelItemCount()):
            item = self.tree.topLevelItem(i)
            if item.text(0) in checked:
                item.setCheckState(0, QtCore.Qt.Checked)
        #get previous results if present
        self.get_results()
        return

    def project_closed(self):
        """Run when parent project is closed"""

        return
