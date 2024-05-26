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
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Phylo
from snipgenie import app, core, widgets, tables, tools, trees
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

    bl = tools.blast_fasta(ref_file, filename, maxseqs=hits)
    bl = get_blast_coverage(bl, ref_file)
    #bl['genus'] = bl.stitle.apply(lambda x: x.split()[1])
    bl['species'] = bl.stitle.apply(lambda x: ' '.join(x.split()[1:3]))
    cols = ['sseqid','sslen','length','perc_cov','pident','stitle','species']
    bl = bl[(bl.pident>=pident) & (bl.perc_cov>=80)].sort_values(['pident','perc_cov'],ascending=False)
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
    menuentry = '16S Species Check'
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

        self.numhits = 20
        self.pident = 95
        self.create_widgets()
        self.fetch_sequence()
        self.results = {}
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

        w = QLabel('Max Hits:')
        vbox.addWidget(w)
        self.hitsentry = w = QSpinBox()
        w.setRange(10,2e3)
        w.setSingleStep(10)
        w.setValue(self.numhits)
        vbox.addWidget(w)
        w = QLabel('Min Percent Identity:')
        vbox.addWidget(w)
        self.pidententry = w = QSpinBox()
        w.setValue(self.pident)
        vbox.addWidget(w)
        w = QLabel('Overwrite:')
        vbox.addWidget(w)
        self.overwritebox = w = QCheckBox()
        w.setChecked(0)
        vbox.addWidget(w)

        button = QPushButton("Show Detailed Results")
        button.clicked.connect(self.show_blast_result)
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
        return

    def run(self):
        """Run against selected genomes"""

        self.parent.opts.applyOptions()
        kwds = self.parent.opts.kwds
        threads = kwds['threads']
        outdir = self.outpath
        hits = self.hitsentry.value()
        pident = self.pidententry.value()
        overwrite = self.overwritebox.isChecked()

        def func(progress_callback):
            table = self.parent.fastq_table
            df = table.model.df.copy()
            rows = table.getSelectedRows()
            new = df.iloc[rows] #subset of sample table to run
            for i,r in new.iterrows():
                if r.assembly is None or not os.path.exists(r.assembly):
                    print ('no assembly found')
                    continue
                sample = r['sample']
                if sample in self.results and overwrite == False:
                    continue
                bl = blast_16S(r.assembly, pident=pident, hits=hits)
                #print (bl)
                self.results[sample] = bl
        self.parent.run_threaded_process(func, self.run_completed)
        return

    def run_completed(self):
        self.parent.processing_completed()
        self.get_results()
        return

    def get_results(self):
        """Reuslts from blast tables"""

        print ('getting results..')
        cols = ['sseqid','perc_cov','pident','stitle','species']
        res=[]
        for sample in self.results:
            bl = self.results[sample]
            x=bl.iloc[0][cols]
            x.name = sample
            res.append(x)

        df = pd.DataFrame(res)
        self.result_table.setDataFrame(df)
        return

    def get_tree(self, bl, name):
        """Get a tree from the top blast results"""

        targetseq = SeqRecord(Seq(bl.iloc[0].qseq), id=name)
        m=bl[:50]
        names = list(m.sseqid)
        #extract relevant sequences from ref fasta
        extract_sequences_by_ids(ref_file, 'temp.fa', names)
        #add target sequence to tree also
        append_sequences_to_fasta('temp.fa',targetseq)
        treefile = get_tree('temp.fa')
        return treefile, m

    def show_blast_result(self):
        """Show detailed results for a sample"""

        idx = self.result_table.getSelectedIndexes()
        sample = idx[0]
        if sample not in self.results:
            return
        cols = ['species','sseqid','perc_cov','pident','length','evalue','bitscore','stitle']
        bl = self.results[sample]
        #print (bl.columns)
        w = QDialog(self.main)
        tabs = QTabWidget(w)
        l = QVBoxLayout()
        w.setLayout(l)
        l.addWidget(tabs)

        t = tables.DataFrameWidget(w, font=core.FONT, fontsize=core.FONTSIZE)
        t.table.setDataFrame(bl[cols])
        idx = tabs.addTab(t, 'Blast Result')

        treefile, m = self.get_tree(bl, sample)
        tv = widgets.TreeViewer(w)
        tv.draw(treefile, df=m.set_index('sseqid'), tiplabelcol='species')
        idx = tabs.addTab(tv, 'Tree')

        sv = widgets.SequencesViewer(w)
        #sv.load_fasta('temp.fa')
        sv.load_alignment('temp.aln')
        #sv.load_records(recs)
        idx = tabs.addTab(sv, 'Sequences')

        w.setGeometry(QtCore.QRect(100, 100, 900, 600))
        w.setWindowTitle(f'{sample} 16S results')
        w.show()
        w.activateWindow()
        return

    def save_data(self):
        """Return save data"""

        data = {}
        data['results'] = self.results
        #data['checked'] = self.get_checked()
        return data

    def load_data(self, data):
        """Load any saved data from project.
        Run when plugin is initially launched."""

        #get previous results if present
        if 'results' in data:
            self.results = data['results']
        self.get_results()
        return

    def project_closed(self):
        """Run when parent project is closed"""

        return
