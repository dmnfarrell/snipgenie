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

import sys,os,subprocess,glob,shutil,re,random
import platform
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Phylo, AlignIO
import numpy as np
import pandas as pd
from  . import tools

qcolors = ['blue','green','crimson','blueviolet','orange','cadetblue','chartreuse','chocolate',
            'coral','gold','cornflowerblue','palegreen','khaki','orange','pink','burlywood',
            'red','lime','mediumvioletred','navy','teal','darkblue','purple','orange',
            'salmon','maroon']

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
    from ete3 import NodeStyle
    for n in t.traverse():
        ns = NodeStyle()
        ns["size"] = 0
        n.set_style(ns)

def remove_tiplabels(t):

    for l in t.iter_leaves():
        l.name = None

def set_nodesize(t, size=12):
    """Change the node size"""

    from ete3 import NodeStyle
    for l in t.iter_leaves():
        clr = l._img_style['fgcolor']
        ns = NodeStyle()
        ns["size"] = size
        #keep color
        ns["fgcolor"] = clr
        l.set_style(ns)

def color_leaves(t, colors, color_bg=False):
    from ete3 import NodeStyle

    for l in t.iter_leaves():
        if l.name in colors:
            #print (l.name, colors[l.name])
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
        if color_bg == True:
            ns["bgcolor"] = clr
        l.set_style(ns)

def get_colormap(values):

    import pylab as plt
    labels = values.unique()
    cmap = plt.cm.get_cmap('Set1')
    colors = [cmap(i) for i in range(len(labels))]
    #colors=qcolors
    #clrs = {labels[i]:cmap(float(i)/(len(labels))) for i in range(len(labels))}
    clrs = dict(list(zip(labels,colors)))
    return clrs

def run_fasttree(infile, bootstraps=100, outpath='.'):
    """Run fasttree"""

    fc = tools.get_cmd('fasttree')
    out = os.path.join(outpath,'fasttree.newick')
    cmd = '{fc} -nt {i} > {o}'.format(fc=fc,b=bootstraps,i=infile,o=out)
    tmp = subprocess.check_output(cmd, shell=True)
    return out

def run_RAXML(infile, name='variants', threads=8, bootstraps=100, outpath='.'):
    """Run Raxml pthreads.
        Returns:
            name of .tree file.
    """

    outpath = os.path.abspath(outpath)
    if not os.path.exists(outpath):
        os.makedirs(outpath, exist_ok=True)

    model = 'GTRCAT'
    s1 = random.randint(0,1e8)
    s2 = random.randint(0,1e8)

    files = glob.glob(os.path.join(outpath,'RAxML_*'))
    for f in files:
        os.remove(f)

    cmd = 'raxmlHPC-PTHREADS -f a -N {nb} -T {t} -m {m} -V -p {s1} -x {s2} -n {n} -w {w} -s {i}'\
            .format(t=threads,nb=bootstraps,n=name,i=infile,s1=s1,s2=s2,m=model,w=outpath)
    print (cmd)
    try:
        tmp = subprocess.check_output(cmd, shell=True)
    except Exception as e:
        print ('Error building tree. Is RAxML installed?')
        return None
    out = os.path.join(outpath,'RAxML_bipartitions.variants')
    return out

def convert_branch_lengths(treefile, outfile, snps):

    tree = Phylo.read(treefile, "newick")
    for parent in tree.find_clades(terminal=False, order="level"):
            for child in parent.clades:
                if child.branch_length:
                    child.branch_length *= snps
    #Phylo.draw(tree)
    Phylo.write(tree, outfile, "newick")
    return

def biopython_draw_tree(filename):

    from Bio import Phylo
    tree = Phylo.read(filename,'newick')
    Phylo.draw(tree)
    return

def create_tree(filename=None, tree=None, ref=None, labelmap=None, colormap=None, color_bg=False, format=1):
    """Draw a tree """

    from ete3 import Tree, PhyloTree, TreeStyle, TextFace
    if filename != None:
        t = Tree(filename, format=format)
    else:
        t = tree
    if ref != None:
        t.set_outgroup(ref)
    if colormap != None:
        color_leaves(t, colormap, color_bg)
    if labelmap != None:
        set_tiplabels(t,labelmap)

    #format_nodes(t)
    ts = TreeStyle()
    return t, ts

def colors_from_labels(df,name,group):
    """Colors from dataframe columns for use with an ete3 tree drawing"""

    labels = df[group].unique()
    colors={}
    i=0
    for l in labels:
        if i>=len(qcolors):
            i=0
        colors[l] = qcolors[i]
        i+=1
    df['color'] = df[group].apply(lambda x: colors[x],1)
    #colormap = dict(zip(df[name],df.color))
    return colors

def remove_nodes(tree, names):

    for n in names:
        node = tree.search_nodes(name=n)[0]
        node.delete()
