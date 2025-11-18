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
import pylab as plt
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
ncbi16S_url = 'https://github.com/dmnfarrell/snipgenie/raw/master/extra/16S_ncbi.fa.gz'
silva_url = 'https://ftp.arb-silva.de/release_138/Exports/SILVA_138_SSURef_NR99_tax_silva.fasta.gz'

def get_blast_coverage(bl, fasta):
    """Get alignment coverage of blast results from original sequence lengths"""

    df=tools.fasta_to_dataframe(fasta)
    df=df.rename(columns={'length':'sslen'})
    bl=bl.merge(df,left_on='sseqid',right_on='name',how='left')
    bl['perc_cov'] = bl.apply(lambda x: round(x.length/x.sslen*100,2),1)
    return bl

def blast_16S_old(filename, database='16S_ncbi.fa', hits=100, pident=99.5):
    """Blast a fasta assembly to 16S sequences"""

    bl = tools.blast_fasta(database, filename, maxseqs=hits)
    bl = get_blast_coverage(bl, database)
    bl = bl.drop_duplicates('sseqid')
    if '16S_ncbi' in database:
        bl['species'] = bl.stitle.apply(lambda x: ' '.join(x.split()[1:3]))
    elif 'SILVA' in database:
        bl['species'] = bl.stitle.apply(lambda x: x.split(';')[-1])
    bl['genus'] = bl.species.apply(lambda x: x.split(' ')[0])
    cols = ['sseqid','sslen','length','perc_cov','pident','stitle','species']
    bl = bl[(bl.pident>=pident) & (bl.perc_cov>=80)].sort_values('pident',ascending=False)
    return bl

def get_blast_coverage_qlen(bl, fasta_query):
    """
    Get alignment coverage of blast results, calculating coverage
    as the alignment length ('length') divided by the query length ('qlen').
    The 'fasta_query' must contain the 16S sequences used as the query.
    """
    # Load the query 16S sequence lengths
    df = tools.fasta_to_dataframe(fasta_query)
    df = df.rename(columns={'length': 'qlen'})
    # Merge the BLAST results with the query sequence lengths
    bl = bl.merge(df, left_on='qseqid', right_on='name', how='left')
    # Calculate percent coverage (alignment length / query length)
    bl['perc_cov'] = bl.apply(lambda x: round(x.length / x.qlen * 100, 2), axis=1)
    return bl

def blast_16S(contig_file, reference, hits=100, pident=99.5):
    """
    Blast 16S reference sequences (query) to contig assembly (target).

    contig_file: The FASTA file of your assembled contigs (the target/subject DB).
    reference: The FASTA file of the 16S reference sequences (the query).
    """

    # --- 1. BLAST: Swap Query and Target ---
    # The 16S reference sequences are the query, the contigs are the database.
    bl = tools.blast_fasta(contig_file, reference, maxseqs=hits)
    # --- 2. Coverage Calculation: Use 16S length (Qlen) ---
    # Calculate coverage based on the 16S reference length (qlen).
    bl = get_blast_coverage_qlen(bl, reference)
    # --- 3. Filtering and Annotation ---
    # The gene is now identified by the query sequence ID (qseqid).
    # We keep the best match for each unique 16S sequence/strain.
    bl = bl.drop_duplicates('qseqid')

    # Parse taxonomic information (assuming 'stitle' still comes from the 16S ref file)
    if '16S_ncbi' in reference:
        # Assuming stitle structure is consistent for NCBI 16S
        #bl['species'] = bl.qseqid.apply(lambda x: ' '.join(x.split()[1:3]))
        bl['species'] = bl.description.apply(lambda x: ' '.join(x.split()[1:3]))
    elif 'SILVA' in reference:
        # Assuming stitle structure is consistent for SILVA
        bl['species'] = bl.description.apply(lambda x: x.split(';')[-1])

    bl['genus'] = bl.species.apply(lambda x: x.split(' ')[0])
    # Filter for high identity and high coverage of the 16S gene
    # perc_cov >= 80 means the alignment covers at least 80% of the 16S reference.
    bl = bl[(bl.pident >= pident) & (bl.perc_cov >= 80)].sort_values('pident', ascending=False)
    # Rename sseqid (contig ID) for clarity in the output
    bl = bl.rename(columns={'sseqid': 'contig_id'})
    cols = ['species', 'genus', 'contig_id', 'qlen', 'length', 'perc_cov', 'pident', 'stitle']
    #print (bl.columns)
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
        self.ref_file = os.path.join(ref_path,'16S_ncbi.fa')
        self.fetch_sequence()
        self.results = {}
        return

    def fetch_sequence(self, url=ncbi16S_url):
        """Get 16S sequences"""

        ref = self.ref_file
        if os.path.exists(self.ref_file):
            return
        print ('fetching 16S sequence data..')
        infile = ref+'.gz'
        import urllib.request
        urllib.request.urlretrieve(url, infile)
        tools.gunzip(infile, ref)
        #tools.make_blast_database(ref)
        return

    def set_database(self):
        """Set current database"""

        name = self.databasew.currentText()
        if name == 'NCBI':
            self.ref_file = os.path.join(ref_path,'16S_ncbi.fa')
            self.fetch_sequence(ncbi16S_url)
        elif name == 'SILVA':
            self.ref_file = os.path.join(ref_path,'SILVA_138_SSURef_NR99_tax_silva.fasta')
            self.fetch_sequence(silva_url)
        return

    def create_widgets(self):
        """Create widgets if GUI 'species', 'genus' plugin"""

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

        w = QLabel('Database:')
        vbox.addWidget(w)
        self.databasew = w = QComboBox()
        w.addItems(['NCBI','SILVA'])
        vbox.addWidget(w)
        w.activated.connect(self.set_database)
        w = QLabel('Max Hits:')
        vbox.addWidget(w)
        self.hitsentry = w = QSpinBox()
        w.setRange(10,2000)
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
            if not 'assembly' in df.columns:
                print('You need to run an assembly for the sample. Use the assembly plugin.')
                return
            rows = table.getSelectedRows()
            new = df.iloc[rows] #subset of sample table to run
            for i,r in new.iterrows():
                sample = r['sample']
                if r.assembly is None or not os.path.exists(r.assembly):
                    print (f'no assembly found for sample {sample}')
                    continue
                if sample in self.results and overwrite == False:
                    continue
                print (f'blasting sample {sample} to 16S..')
                tools.make_blast_database(r.assembly)
                bl = blast_16S(r.assembly, self.ref_file, pident=pident, hits=hits)
                #remove temp files

                #print (bl)
                if len(bl) == 0:
                    print ('no blast hits found!')
                else:
                    self.results[sample] = bl
        self.parent.run_threaded_process(func, self.run_completed)
        return

    def run_completed(self):
        self.parent.processing_completed()
        self.get_results()
        return

    def get_results(self):
        """Results from blast tables"""

        print ('getting results..')
        cols = ['species','perc_cov','pident','stitle']#,'genus']
        res=[]
        for sample in self.results:
            bl = self.results[sample]
            x = bl.iloc[0][cols]
            x.name = sample
            res.append(x)
            #print (res)
        df = pd.DataFrame(res)
        self.result_table.setDataFrame(df)
        return

    def get_tree(self, bl, name):
        """Get a tree from the top blast results"""

        targetseq = SeqRecord(Seq(bl.iloc[0].qseq), id=name)
        m=bl[:50]
        names = list(m.qseqid)
        #extract relevant sequences from ref fasta
        extract_sequences_by_ids(self.ref_file, 'temp.fa', names)
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
        cols = ['species','qseqid','perc_cov','pident','length','evalue','bitscore','stitle']
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
        tv.draw(treefile, df=m.set_index('qseqid'), tiplabelcol='species')
        idx = tabs.addTab(tv, 'Tree')

        sv = widgets.SequencesViewer(w)
        #sv.load_fasta('temp.fa')
        sv.load_alignment('temp.aln')
        #sv.load_records(recs)
        idx = tabs.addTab(sv, 'Sequences')

        #fig,ax = plt.subplots(1,1)
        #pw = widgets.PlotWidget(figure=fig)
        #bl.genus.value_counts().plot(kind='barh',ax=ax)
        #plt.tight_layout()
        #idx = tabs.addTab(pw, 'Genus Distribution')

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
