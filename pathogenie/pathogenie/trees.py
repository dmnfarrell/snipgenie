"""
    Tree methods for bacterial phylogenetics, mostly using ete3.
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

from __future__ import print_function
import sys,os,subprocess,glob,shutil,re
import platform
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Phylo, AlignIO
import numpy as np
import pandas as pd
sys.path.append('ete')
from ete3 import Tree, NodeStyle, TreeStyle

def set_tiplabels(t, labelmap):
    for l in t.iter_leaves():
        #print (l.name)
        if l.name in labelmap:            
            l.name = labelmap[l.name]
    return    

def delete_nodes(t, names):
    for l in t.iter_leaves():        
        if l.name in names:          
            l.delete()   
            
def format_nodes(t):
    for n in t.traverse(): 
        ns = NodeStyle()
        ns["size"] = 0
        n.set_style(ns)
        
def color_leaves(t, colors):
    for l in t.iter_leaves():  
        if l.name in colors:
            clr = colors[l.name]
        else:
            clr='black'
        #print (clr)
        # create a new label with a color attribute
        #N = AttrFace("name", fgcolor=clr)
        #l.add_face(N, 1, position='aligned')
        ns = NodeStyle()
        ns["size"] = 12
        ns["fgcolor"] = clr
        l.set_style(ns)
        
def get_colormap(values):
    labels = values.unique()
    cmap = plt.cm.get_cmap('Set1')    
    colors=qcolors
    #clrs = {labels[i]:cmap(float(i)/(len(labels))) for i in range(len(labels))}
    clrs = dict(list(zip(labels,colors)))
    return clrs   