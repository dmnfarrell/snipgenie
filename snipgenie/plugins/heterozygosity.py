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

class HeteroCheckerPlugin(Plugin):
    """Hetero checker plugin for SNiPgenie"""

    #uncomment capabilities list to appear in menu
    capabilities = ['gui','docked']
    requires = ['']
    menuentry = 'Heterozygosity Check'
    name = 'Heterozygosity Check'
    iconfile = 'heterozygosity.svg'
    side = 'right' #dock location to load plugin

    def __init__(self, parent=None):
        """Customise this and/or create_widgets"""

        if parent==None:
            return
        self.parent = parent
        #self.outpath = os.path.join(self.parent.outputdir, 'contam')
        #if not os.path.exists(self.outpath):
        #    os.makedirs(self.outpath, exist_ok=True)
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

        #t = self.result_table = tables.DataFrameTable(self.main)
        t = self.table_widget = tables.DataFrameWidget(self.main, toolbar=True, app=self.parent)
        self.result_table = self.table_widget.table
        layout.addWidget(t)
        bw = self.create_buttons(self.main)
        layout.addWidget(bw)
        return

    def create_buttons(self, parent):

        bw = QWidget(parent)
        vbox = QVBoxLayout(bw)

        w = QLabel('Pivot result:')
        vbox.addWidget(w)
        self.pivotbox = w = QCheckBox()
        w.setChecked(1)
        vbox.addWidget(w)
        w = QLabel('Value column:')
        vbox.addWidget(w)
        self.valuecolw = w = QComboBox()
        w.addItems(['het','DP','AD','ADF','ADR','QUAL','ALT','mut'])
        vbox.addWidget(w)
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

    def load_table(self):
        """load genomes table"""

        self.refs = pd.read_csv(self.ref_table, index_col=0)
        self.update_treelist()
        return

    def run(self):
        """Run against selected genomes"""

        def het(x):
            if sum(x.AD) == 0:
                return
            return round(min(x.AD)/sum(x.AD),3)

        result = []
        def func(progress_callback):
            table = self.parent.fastq_table
            df = table.model.df.copy()
            rows = table.getSelectedRows()
            data = df.iloc[rows]
            samples = list(data['sample'].unique())
            #old way
            vcf_file = os.path.join(self.parent.outputdir, 'merged.vcf.gz')
            vdf = tools.vcf_to_dataframe(vcf_file)

            i=0
            for s in samples:
                print (s)
                #vcf_file = os.path.join(self.parent.outputdir, 'variant_calling', s, 'snps.bcf')
                #vdf = tools.vcf_to_dataframe(vcf_file)
                x = vdf[vdf['sample']==s].copy()
                x['het'] = x.apply(het,1)
                i+=1
                h = x[x['het']>0.1]
                #sites.append(h)
                result.append(h)
            #print (sites)
            self.results = pd.concat(result)

        self.parent.run_threaded_process(func, self.run_completed)
        return

    def run_completed(self):
        self.parent.processing_completed()
        self.get_results()
        return

    def get_results(self):
        """Reuslts from mapping"""

        print ('getting results..')
        df = self.results
        #print (df[:10])
        if self.pivotbox.isChecked():
            valuecol = self.valuecolw.currentText()
            p = pd.pivot_table(df, index='pos',columns='sample',values=valuecol,aggfunc='first')
        else:
            p = df
        self.result_table.setDataFrame(p)
        return

    def plot_results(self):
        return

    def save_data(self):
        """Return save data"""

        data = {}
        data['result'] = self.result_table.model.df
        return data

    def load_data(self, data):
        """Load any saved data from project.
        Run when plugin is initially launched."""

        if 'result' in data:
            self.result_table.setDataFrame(data['result'])
        return

    def project_closed(self):
        """Run when parent project is closed"""

        return
