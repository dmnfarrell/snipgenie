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

def snp_table_from_matrix(nucmat):
    """Re-format core snp matrix into long form dataframe"""

    nucmat = nucmat.set_index('pos').T.stack().reset_index()
    nucmat.columns=['sample','pos','nuc']
    nucmat['allele'] = nucmat.pos.astype(str)+nucmat.nuc
    return nucmat

def get_snps(nucmat, sample):
    """Get sample snps from a core snp matrix of samples.
    This matrix is the output from snipgenie (core.txt)"""

    df = get_snp_table(nucmat)
    return list(df[df['sample']==sample].allele)

def get_clade_snps(nucmat, clusts, col='snp12'):
    """Get clade specific snps from the snp matrix.
    Args:
        nucmat: core snp matrix from snipgenie
        clusts: dataframe of sequence clusters with varying snp cutoffs.
            generated using trees.get_clusters(tree file)
        col: snp level
    """

    c = snp_table_from_matrix(nucmat)
    snp12cl = clusts[~clusts.snp12.isin([-1,1])]
    c = c.merge(snp12cl,left_on='sample',right_on='SequenceName')
    g = c.groupby(['pos',col,'allele']).agg({'allele':np.size})
    g = g.rename(columns={'allele':'size'}).reset_index()
    #print (g)
    found=[]
    for i,df in g.groupby('pos'):
        vc = df.allele.value_counts()
        #get only alleles presents in one cluster
        vc = (vc[vc==1])
        if len(vc)==1:
            a = vc.index[0]
            f = df[df.allele==a]
            found.append(df[df.allele==a])
    found=pd.concat(found)
    #print (found[found.snp12==10])
    found.groupby(col).size()
    return found

def lookup_sample(clade_snps, x):
    """Lookup a samples snps to identify known clade"""

    found = clade_snps[clade_snps.allele.isin(x)]
    #print (found)
    found = set(found.snp12)
    if len(found) == 1:
        return list(found)[0]
    elif len(found)==0:
        return
    return found

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
