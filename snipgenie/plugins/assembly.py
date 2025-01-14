"""
    Genome assembly plugin for snipgenie
    Created January 2024
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

def calculate_N50(list_of_lengths):
    """Calculate N50 for a sequence of numbers.
    Args:
        list_of_lengths (list): List of numbers.
    Returns:
        float: N50 value.
    """
    tmp = []
    for tmp_number in set(list_of_lengths):
            tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
    tmp.sort()
    if (len(tmp) % 2) == 0:
        median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
    else:
        median = tmp[int(len(tmp) / 2)]
    return median

def get_sequence_lengths(fasta_file):
    """Get sequence lengths from fasta"""

    lengths = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        lengths.append(len(record.seq))
    return lengths

class AssemblyPlugin(Plugin):
    """Genome assembly plugin for SNiPgenie"""

    #uncomment capabilities list to appear in menu
    capabilities = []
    requires = []
    menuentry = 'Genome Assembly'
    name = 'Genome Assembly'
    iconfile = 'assembly.svg'

    def __init__(self, parent=None):
        """Customise this"""

        if parent==None:
            return
        self.parent = parent
        return

    def run(self):
        """Run assembly, updates parent table in place."""

        table = self.parent.fastq_table
        df = table.model.df
        rows = table.getSelectedRows()
        data = df.iloc[rows]
        self.parent.opts.applyOptions()
        kwds = self.parent.opts.kwds
        threads = kwds['threads']
        path = os.path.join(self.parent.outputdir, 'assembly')
        #path for final files
        fastapath = os.path.join(self.parent.outputdir, 'assembled')
        if not os.path.exists(fastapath):
            os.makedirs(fastapath,exist_ok=True)

        if data is None or len(data) == 0:
            return
        msg = 'This will de novo assemble the reads.\nIt may take some time. Are you sure?\n'
        reply = QMessageBox.question(self.parent, 'Warning!', msg,
                                            QMessageBox.No | QMessageBox.Yes )
        if reply == QMessageBox.No:
            return

        print ('Running de novo assembly. This may take some time.')
        def func(progress_callback):
            res=[]
            cols = ['sample','spotype','sb']
            for i,r in data.iterrows():
                name = r['sample']
                outfile = os.path.join(fastapath,'%s.fa') %name
                if not os.path.exists(outfile):
                    tools.spades(r.filename1,r.filename2, os.path.join(path,name),
                                  outfile, threads)
                else:
                    print (f'assembly exists in {outfile}, remove this file to recreate')
                df.loc[i,'assembly'] = outfile
                lengths = get_sequence_lengths(outfile)
                df.loc[i,'N50'] = calculate_N50(lengths)
            return

        self.parent.run_threaded_process(func, self.parent.processing_completed)
        return

    def run_completed(self):
        self.parent.processing_completed()
        return
