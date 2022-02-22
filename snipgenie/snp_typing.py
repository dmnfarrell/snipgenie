"""
    SNP typing module for snipgenie.
    Created Mar 2021
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

import sys,os,subprocess,glob,shutil,re,random,time
import numpy as np
import pandas as pd
import pylab as plt
import scipy.cluster.hierarchy as shc
from sklearn.preprocessing import normalize
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO,AlignIO
from Bio import Phylo
import seaborn as sns
import toytree
from snipgenie import app, trees, tools

home = os.path.expanduser("~")
config_path = os.path.join(home,'.config','snipgenie')
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'data')

#default bovis snp matrix
nucmat = pd.read_csv(os.path.join(datadir,'nuc_snps_ireland.txt'),sep=' ')
clusts = pd.read_csv(os.path.join(datadir,'ireland_clusters.txt'),sep='\t')
clade_snps = pd.read_csv(os.path.join(datadir,'ireland_clade_snps.csv'))

def run_tree_cluster(f,dist):

    cmd = 'TreeCluster.py  -i {f} -t {d}'.format(f=f,d=dist)
    cl=subprocess.check_output(cmd, shell=True)
    cl=pd.read_csv(io.BytesIO(cl),sep='\t')
    return cl

def snps_to_fasta(snpmat, outfile):
    """Write snp matrix to fasta file"""

    snpmat = snpmat.fillna('N')
    recs = []
    for col in snpmat.columns:
        seq = ''.join(snpmat[col])
        seqrec = SeqRecord(Seq(seq),id=col)
        recs.append(seqrec)
    SeqIO.write(recs, outfile, 'fasta')
    return

def tree_from_snps(snpmat):
    """Make tree from snp matrix"""

    snps_to_fasta(snpmat, 'snps.fa')
    treefile = trees.run_fasttree('snps.fa')
    tre = toytree.tree(treefile)
    mystyle = { "layout": 'r','node_sizes':1,'tip_labels_align':False}
    tre.ladderize().draw(**mystyle,width=700)
    return tre

def make_ref_snps(nucmat, clusts, column='ClusterNumber'):
    """Add cluster info to snps"""

    nucmat = nucmat.set_index('pos')
    X=nucmat.T.merge(clusts,left_index=True,right_on='SequenceName').set_index([column]).T
    return X

def get_clade_snps(refmat):
    """Get unique clade SNPs from a SNP matrix
       returns: a dataframe with unique positions/allele for each clade
       with this format
              clade      pos allele
           2   490878      G
           2   804997      T
           2   941068      A
           2  1124266      G
    """

    res=[]
    clusters = refmat.columns.unique()
    for c in clusters:
        for pos,r in list(refmat.iterrows())[:700]:
            #print (pos)
            a = r[c]
            b = r[~r.index.isin([c])]
            f1 = a.value_counts()
            f2 = b.value_counts()
            alt1 = f1.index[0]
            if len(f1)>1:
                continue
            alt2 = f2.index[0]
            if alt1 in f2:
                continue
            res.append((c,pos,alt1))

    res = pd.DataFrame(res,columns=['clade','pos','allele'])
    print (res)
    return res

def lookup_sample(snptable, snps):
    """Look up a sample using snps and known clades
        snptable: reference lookup table
        snps: a series with snps at each position for the
        given sample, this can be derived from a single row
        in the snp matrix produced from snipgenie
    """

    found=[]
    for i,r in snptable.iterrows():
        if not r.pos in snps.index:
            continue
        if snps[r.pos] == r.allele:
            #print (r.pos,r.allele,r.clade)
            found.append(r.clade)
    if len(found) == 0:
        return
    return set(found)

def type_samples(nucmat):
    """
    Type multiple samples.
    Args:
        nucmat: a dataframe with the following format-
        pos       687  937  1303 ..
        sample1    C    A    G
        sample2    C    A    G
        ...
    Returns:
        types for each sample
    """
    snptable = clade_snps
    for name,r in nucmat.iterrows():
        #print (r)
        cl = lookup_sample(snptable, r)
        print (name,cl)

def encode_snps(x):
    """encode snps as string for storage"""

    s=[]
    for i in zip(x.index.astype(str),x.values):
        s.append(''.join(i))
    s = ';'.join(s)
    return s

def decode_snps(s):
    """decode snps"""

    x=s.split(';')
    pos=[]
    alleles=[]
    for i in x:
        n,p,a = re.split(r'(\d+)', i)
        pos.append(p)
        alleles.append(a)
    x = pd.Series(alleles,pos)
    x.index.name='pos'
    return x

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='snipgenie typing tool.')
    parser.add_argument("-i", "--input", action='append', dest="input", default=[],
                        help="input ", metavar="FILE")


if __name__ == '__main__':
    main()
