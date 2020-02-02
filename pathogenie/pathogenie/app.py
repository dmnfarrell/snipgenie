#!/usr/bin/env python

"""
    pathogenie cmd line tool.
    Created Nov 2019
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 3
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

from __future__ import absolute_import, print_function
import sys,os,subprocess,glob,re
import urllib, hashlib, shutil
import tempfile
import pandas as pd
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna
from . import tools

tempdir = tempfile.gettempdir()
home = os.path.expanduser("~")
config_path = os.path.join(home,'.config/pygenefinder')
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'data')
dbdir = os.path.join(config_path, 'db')
customdbdir = os.path.join(config_path, 'custom')

def assign_sample_ids(filenames):
    """Assign new ids for sample filenames e.g. files. Useful for
       replacing long file names with short ids.
       Returns: dataframe with labels
    """

    i=1
    labels = {}
    for fname in filenames:
        n = os.path.splitext(os.path.basename(fname))[0]
        print (n)
        sid = 's%02d' %i
        labels[n] = sid
        i+=1
    l = pd.DataFrame.from_dict(labels,orient='index')
    l.columns = ['id']; l.index.name='filename'
    return l
