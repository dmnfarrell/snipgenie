"""
    Spotpying plugin for snipgenie
    Created November 2022
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

class SpoligotypingPlugin(Plugin):
    """Contam checker plugin for SNiPgenie"""

    #uncomment capabilities list to appear in menu
    capabilities = []
    requires = []
    menuentry = 'M.bovis Spoligotyping'
    name = 'M.bovis Spoligotyping'
    iconfile = 'bar-code.svg'

    def __init__(self, parent=None):
        """Customise this"""

        if parent==None:
            return
        self.parent = parent
        return

    def run(self):
        """Run spoligotyping, updates parent table in place."""

        print ('running spoligotyping..')
        table = self.parent.fastq_table
        df = table.model.df
        rows = table.getSelectedRows()
        data = df.iloc[rows]
        self.parent.opts.applyOptions()
        kwds = self.parent.opts.kwds

        if data is None or len(data) == 0:
            return

        def func(progress_callback):
            res=[]
            cols = ['sample','spotype','sb']
            for i,r in data.iterrows():
                name = r['sample']
                s = tools.get_spoligotype(r.filename1, threads=kwds['threads'])
                sb = tools.get_sb_number(s)
                print (name, s, sb)
                #set new values in place
                df.loc[i,cols] = [name,s,sb]
            return

        self.parent.run_threaded_process(func, self.run_completed)
        return

    def run_completed(self):
        self.parent.processing_completed()
        return
