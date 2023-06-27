"""
    Strain typing plugin for snipgenie
    Created May 2023
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
import pylab as plt
from Bio import Phylo, SeqIO, Entrez
import toytree, toyplot
from snipgenie import app, widgets, tables, tools, trees, clustering, plotting, treeview
from snipgenie.plugin import Plugin

#location for main snp tables/clustering

datadir = os.path.join(app.module_path,'mbovis_typing')
core_snps = os.path.join(datadir, 'core_snps_mbovis.txt')
cluster_members = os.path.join(datadir, 'cluster_members.parquet')
core_snp_dist = os.path.join(datadir, 'snpdist.csv')

class ResultsViewer(QWidget):
    """Widget for results plots"""
    def __init__(self, parent=None, meta=None):

        super(ResultsViewer, self).__init__(parent)
        self.setGeometry(QtCore.QRect(200, 200, 900, 600))
        self.main = QTabWidget()
        self.setLayout(QVBoxLayout())
        self.layout().addWidget(self.main)
        self.mstplot = widgets.PlotWidget()
        self.main.addTab(self.mstplot,'mst')
        self.tv = treeview.TreeViewer(meta=meta)
        self.main.addTab(self.tv,'phylogeny')
        self.clustermap = widgets.PlotWidget()
        self.main.addTab(self.clustermap,'clustermap')
        #self.distmatrix = tables.DistMatrixTable(fontsize=12)
        self.distmatrix = tables.DataFrameTable()
        self.main.addTab(self.distmatrix,'dist matrix')
        self.snptable = tables.SNPTable()
        self.main.addTab(self.snptable,'SNPs')
        self.setWindowTitle('strain clustering')
        return

class StrainTypingPlugin(Plugin):
    """M.bovis strain typing plugin for SNiPgenie"""

    #uncomment capabilities list to appear in menu
    capabilities = ['gui','docked']
    requires = ['']
    menuentry = 'M.bovis Strain Typing'
    name = 'M.bovis Typing'
    iconfile = 'mbovis.svg'
    side = 'right'

    def __init__(self, parent=None):
        """Customise this and/or create_widgets"""

        if parent==None:
            return
        self.parent = parent
        self.outpath = os.path.join(self.parent.outputdir, 'mbovis_typing')
        os.makedirs(self.outpath, exist_ok=True)
        self.create_widgets()
        #load previous
        stfile = os.path.join(self.outpath,'strains.csv')
        snpdistfile = os.path.join(self.outpath,'new_snpdist.csv')
        snpfile = os.path.join(self.outpath,'new_snps.csv')
        if os.path.exists(stfile):
            try:
                self.st = pd.read_csv(stfile, index_col=0)
                self.result_table.setDataFrame(self.st)
                self.newdm = pd.read_csv(snpdistfile, index_col=0)
                self.newsnps = pd.read_csv(snpfile, index_col=0, sep=' ')
                self.result_table.setDataFrame(self.st)
            except:
                pass

        self.meta = self.parent.meta
        return

    def create_widgets(self):
        """Create widgets if GUI plugin"""

        self.main = QWidget()
        layout = self.layout = QVBoxLayout()
        layout.setAlignment(QtCore.Qt.AlignTop)
        self.main.setLayout(layout)

        t = self.table_widget = tables.DataFrameWidget(self.main, app=self.parent)
        self.result_table = self.table_widget.table
        layout.addWidget(t)
        bw = self.create_buttons(self.main)
        layout.addWidget(bw)
        return

    def create_buttons(self, parent):
        """Add buttons"""

        bw = QWidget(parent)
        vbox = QVBoxLayout(bw)

        button = QPushButton("Reset")
        button.clicked.connect(self.show_all)
        vbox.addWidget(button)

        row = QWidget()
        vbox.addWidget(row)
        hb = QHBoxLayout(row)
        button = QPushButton("Find closest samples")
        button.clicked.connect(lambda : self.find_samples('closest'))
        hb.addWidget(button)
        w = QLabel('No. Samples:')
        hb.addWidget(w)
        self.closestentry = w = QLineEdit()
        w.setText('8')
        w.resize(20, 20)
        hb.addWidget(w)

        row = QWidget()
        vbox.addWidget(row)
        hb = QHBoxLayout(row)
        button = QPushButton("Find within n SNPs")
        button.clicked.connect(lambda : self.find_samples('within'))
        hb.addWidget(button)
        w = QLabel('threshold:')
        hb.addWidget(w)
        self.snpthresholdentry = w = QLineEdit()
        w.setText('3')
        w.resize(20, 20)
        hb.addWidget(w)

        button = QPushButton("Plot results")
        button.clicked.connect(self.plot_results)
        vbox.addWidget(button)
        button = QPushButton("Run typing")
        button.clicked.connect(self.run)
        vbox.addWidget(button)
        return bw

    def show_selected(self):
        """Show only selected rows"""

        r = self.result_table.getSelectedRows()
        samples = self.st.iloc[r]
        self.idx = df.index
        self.result_table.setDataFrame(samples)
        return

    def find_samples(self, how='closest'):
        """Show closest related members of selected sample"""

        r = self.result_table.getSelectedRows()
        if len(r) == 0:
            return
        df = self.st.iloc[r]
        s = df.index[0]

        if how == 'closest':
            n = int(self.closestentry.text())
            selected = tools.get_closest_samples(self.newdm, s, n)
        elif how == 'within':
            n = int(self.snpthresholdentry.text())
            selected = tools.get_within_distance(self.newdm, s, n)
        if len(selected) <= 1:
            self.idx = None
            print ('no samples found')
            return
        self.idx = list(selected.index)+[s]
        samples = self.st.loc[self.idx]
        self.result_table.setDataFrame(samples)
        return

    def plot_results(self):
        """Plot current selection"""

        idx = self.idx
        dm = self.newdm.loc[idx,idx]
        samples = self.st.loc[idx]
        self.result_table.setDataFrame(samples)

        if not hasattr(self, 'resultsview'):
            self.resultsview = ResultsViewer(meta=self.meta)

        #phylogeny
        snpmat = self.newsnps[['pos']+idx]
        #print(snpmat)
        treefile = trees.tree_from_snps(snpmat)
        self.resultsview.tv.load_tree(treefile)
        self.resultsview.tv.update()

        #min spanning tree
        self.resultsview.show()
        self.resultsview.activateWindow()
        ax = self.resultsview.mstplot.ax
        tools.dist_matrix_to_mst(dm, ax)
        #self.plotview.canvas.draw()

        #dist matrix
        self.resultsview.distmatrix.setDataFrame(dm)

        #snps
        c = self.newsnps.set_index('pos').T
        c = c.loc[idx]
        c = c[[i for i in c if c[i].nunique()>1]]
        #print(c)
        self.resultsview.snptable.setDataFrame(c)

        #clustermap
        #import seaborn as sns
        #g = sns.clustermap(dm,xticklabels=True,yticklabels=True,annot=True,fmt='g',cmap='Blues')
        #plt.tight_layout()
        #self.resultsview.clustermap.figure = g.figure
        return

    def show_all(self):
        """show all strains"""

        df = self.st
        self.result_table.setDataFrame(df)
        return

    def run_typing(self):
        """Run typing of samples in project"""

        stfile = os.path.join(self.outpath,'strains.csv')
        if os.path.exists(stfile):
            answer = QMessageBox.question(self, 'Continue?',
                             'Overwrite last results?',
                             QMessageBox.Yes, QMessageBox.No)
            if not answer:
                return
        def func(progress_callback):

            #rows = table.getSelectedRows()
            #new = df.iloc[rows] #subset of sample table to run

            coresnps = pd.read_csv(core_snps, sep=' ')
            members = pd.read_parquet(cluster_members)
            snpdist = pd.read_csv(core_snp_dist,index_col=0)
            print('core snps loaded: %s x %s samples' %(len(coresnps),len(coresnps.columns)))

            #get new snps from current project
            currsnps = pd.read_csv(self.parent.snp_matrix, sep=' ',index_col=0)
            current = coresnps.index

            print ('combining new snps..')
            self.newsnps = df = tools.combine_core_snps(coresnps, currsnps)
            print ('getting new alignment..')
            aln = tools.alignment_from_snps(df)
            print (len(aln))
            print ('computing new distance matrix..')
            newdm = tools.update_snp_dist_matrix(aln, snpdist)

            clusts,newmembers = clustering.get_cluster_levels(newdm, members)
            newst = clustering.generate_strain_names(clusts)
            print ('done')
            #print (newst[:10])
            self.st = newst
            self.newdm = newdm

        self.parent.run_threaded_process(func, self.run_completed)
        return

    def run_completed(self):
        self.parent.processing_completed()
        self.get_results()
        return

    def get_results(self):

        #table = self.parent.fastq_table
        #df = table.model.df.copy()
        self.result_table.setDataFrame(self.st)
        #save current table
        self.st.to_csv(os.path.join(self.outpath,'strains.csv'))
        self.newdm.to_csv(os.path.join(self.outpath,'new_snpdist.csv'))
        self.newsnps.to_csv(os.path.join(self.outpath,'new_snps.csv'), sep=' ')
        return

    def project_closed(self):
        """Run when parent project is closed"""

        return
