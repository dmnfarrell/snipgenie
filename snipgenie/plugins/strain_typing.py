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
import pylab as plt
from Bio import Phylo, SeqIO
from Bio import Entrez
from snipgenie import app, widgets, tables, tools, clustering, plotting
from snipgenie.plugin import Plugin

#location for main snp tables/clustering

datadir = os.path.join(app.module_path,'strain_typing')
core_snps = os.path.join(datadir, 'core_snps_mbovis.txt')
cluster_members = os.path.join(datadir, 'cluster_members.parquet')
core_snp_dist = os.path.join(datadir, 'snpdist.csv')

class StrainTypingPlugin(Plugin):
    """M.bovis strain typing plugin for SNiPgenie"""

    #uncomment capabilities list to appear in menu
    capabilities = ['gui','docked']
    requires = ['']
    menuentry = 'M.bovis Strain Typing'
    name = 'Strain Typing'
    iconfile = 'mbovis.svg'
    side = 'right' 

    def __init__(self, parent=None):
        """Customise this and/or create_widgets"""

        if parent==None:
            return
        self.parent = parent  
        self.outpath = os.path.join(self.parent.outputdir, 'strain_typing')
        os.makedirs(self.outpath, exist_ok=True)
        self.create_widgets()
        #load previous
        stfile = os.path.join(self.outpath,'strains.csv')
        snpdistfile = os.path.join(self.outpath,'new_snpdist.csv')
        if os.path.exists(stfile):
            self.st = pd.read_csv(stfile, index_col=0)
            self.result_table.setDataFrame(self.st)
            self.newdm = pd.read_csv(snpdistfile, index_col=0)
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

        bw = QWidget(parent)
        vbox = QVBoxLayout(bw)
        button = QPushButton("Run current")
        button.clicked.connect(self.run)
        vbox.addWidget(button)
        button = QPushButton("Show closest samples")
        button.clicked.connect(self.show_closest)
        vbox.addWidget(button)
        w = QLabel('Samples to find:')
        vbox.addWidget(w)
        self.closestentry = w = QLineEdit()
        w.setText('5')
        vbox.addWidget(w)     
        return bw

    def show_closest(self):
        """Show related members of selected sample"""

        r = self.result_table.getSelectedRows()
        df = self.st.iloc[r]
        s = df.index[0]
        #print (s)     
        n = int(self.closestentry.text())
        close = tools.get_closest_samples(self.newdm, s, n)
        print (close)
        x = list(close.index)+[s]
        dm = self.newdm.loc[x,x]
        print (dm)
        #treefile = trees.run_fasttree(outfasta, self.outdir)

        #min spanning tree        
        self.plotview = widgets.PlotWidget()
        self.plotview.show()
        #self.plotview.activateWindow()  
        ax = self.plotview.ax
        tools.dist_matrix_to_mst(dm, ax)
        #import seaborn as sns
        #g=sns.clustermap(dm,xticklabels=True,yticklabels=True,annot=True,cmap='Blues')
        #plt.tight_layout()
        #self.plotview.figure = g.figure
        return
    
    def cluster_map(self):

        return
    
    def run(self):
        """Run typing of samples in project"""
        
        def func(progress_callback):

            #rows = table.getSelectedRows()
            #new = df.iloc[rows] #subset of sample table to run
   
            coresnps = pd.read_csv(core_snps, sep=' ')
            members = pd.read_parquet(cluster_members)
            snpdist = pd.read_csv(core_snp_dist,index_col=0)
            print('core snps loaded: %s x %s samples' %(len(coresnps),len(coresnps.columns)))
           
            #get new snps from current project
            newsnps = pd.read_csv(self.parent.snp_matrix, sep=' ',index_col=0)
            current = coresnps.index            

            print ('combining new snps..')
            df = tools.combine_core_snps(coresnps, newsnps)
            print ('getting new alignment..')
            aln = tools.alignment_from_snps(df)
            print (len(aln))
            print ('computing new distance matrix..')
            newdm = tools.update_snp_dist_matrix(aln, snpdist)

            clusts,newmembers = clustering.get_cluster_levels(newdm, members)
            newst = clustering.generate_strain_names(clusts)
            print 
            print (newst[:10])
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
        
        return

    def project_closed(self):
        """Run when parent project is closed"""

        return
