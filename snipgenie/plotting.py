"""
    Plotting methods for snpgenie
    Created Jan 2020
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
import os, sys, io, random, subprocess, time
import string
import numpy as np
import pandas as pd
pd.set_option('display.width', 200)
from importlib import reload

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO
from pyfaidx import Fasta
import pylab as plt
from . import tools

def heatmap(df, cmap='gist_gray_r', w=15, h=5, ax=None):
    """Plot dataframe matrix"""

    if ax == None:
        fig, ax = plt.subplots(figsize=(w,h))
    im = ax.pcolor(df, cmap=cmap)
    ax.set_xticks(np.arange(len(df.columns))+0.5)
    ax.set_yticks(np.arange(len(df))+0.5)
    ax.set_xticklabels(df.columns)
    ax.set_yticklabels(df.index)
    plt.setp(ax.get_xticklabels(), rotation=80, ha="right",
             rotation_mode="anchor")
    plt.tight_layout()
    return

def get_fasta_length(filename, key=None):
    """Get length of reference sequence"""

    refseq = Fasta(filename)
    if key==None:
        key = list(refseq.keys())[0]
    l = len(refseq[key])
    return l

def get_fasta_names(filename):
    """Get names of fasta sequences"""

    refseq = Fasta(filename)
    return list(refseq.keys())

def get_fasta_sequence(filename, start, end, key=0):
    """Get chunk of indexed fasta sequence at start/end points"""

    from pyfaidx import Fasta
    refseq = Fasta(filename)
    if type(key) is int:
        chrom = list(refseq.keys())[key]
    print (chrom)
    seq = refseq[chrom][start:end].seq
    return seq

def get_chrom(bam_file):
    """Get first sequence name in a bam file"""

    import pysam
    samfile = pysam.AlignmentFile(bam_file, "r")
    iter=samfile.fetch(start=0,end=10)
    for read in iter:
        if read.reference_name:
            return read.reference_name

def get_coverage(bam_file, chr, start, end):
    """Get coverage from bam file at specified region"""

    import pysam
    if bam_file is None or not os.path.exists(bam_file):
        return
    samfile = pysam.AlignmentFile(bam_file, "r")
    vals = [(pileupcolumn.pos, pileupcolumn.n) for pileupcolumn in samfile.pileup(chr, start, end)]
    df = pd.DataFrame(vals,columns=['pos','coverage'])
    df = df[(df.pos>=start) & (df.pos<=end)]
    #fill with zeroes if there is no data at ends
    if df.pos.max() < end:
        new = pd.DataFrame({'pos':range(df.pos.max(), end)})
        new['coverage'] = 0
        df = df.append(new).reset_index(drop=True)
    return df

def get_bam_aln(bam_file, chr, start, end, group=False):
    """Get all aligned reads from a sorted bam file for within the given coords"""

    import pysam
    if not os.path.exists(bam_file):
        return
    if chr is None:
        return
    if start<1:
        start=0
    samfile = pysam.AlignmentFile(bam_file, "r")
    iter = samfile.fetch(chr, start, end)
    d=[]
    for read in iter:
        st = read.reference_start
        d.append([read.reference_start, read.reference_end, read.cigarstring,
                  read.query_name,read.query_length,read.mapping_quality])
    df = pd.DataFrame(d,columns=['start','end','cigar','name','length','mapq'])
    if len(df) == 0:
        return pd.DataFrame()
    if group == True:
        df['counts'] = df.groupby(['start','end']).name.transform('count')
        df = df.drop_duplicates(['start','end'])
    df['y'] = 1
    bins = (end-start)/150
    if bins < 1:
        bins = 1
    xbins = pd.cut(df.start,bins=bins)
    df['y'] = df.groupby(xbins)['y'].transform(lambda x: x.cumsum())
    return df

def plot_coverage(df, plot_width=800, plot_height=60, xaxis=True, ax=None):
    """Plot a bam coverage dataframe returned from get_coverage
    Args:
        df: dataframe of coverage data (from get_coverage)
        plot_width: width of plot
        xaxis: plot the x-axis ticks and labels
    """

    #if df is None or len(df)==0:
    #    return plot_empty(plot_width=plot_width,plot_height=plot_height)
    df['y'] = df.coverage/2
    x_range = (df.pos.min(),df.pos.max())
    top = df.coverage.max()
    if ax==None:
        fig,ax = plt.subplots(1,1,figsize=(15,1))
    ax.fill_between(df.pos,df.y,color='gray')
    ax.set_xlim(x_range)
    if xaxis==False:
        ax.get_xaxis().set_visible(False)
    return

def plot_bam_alignment(bam_file, chr, xstart, xend, ystart=0, yend=100,
                        rect_height=.6, fill_color='gray', ax=None):
    """bam alignments plotter.
    Args:
        bam_file: name of a sorted bam file
        start: start of range to show
        end: end of range
    """

    h = rect_height
    #cover the visible range from start-end
    o = (xend-xstart)/2
    #get reads in range into a dataframe
    df = get_bam_aln(bam_file, chr, xstart-o, xend+o)
    #print (df[:4])
    df['x'] = df.start+df.length/2
    df['y'] = df.y*(h+1)
    #set colors by quality
    df['color'] = df.apply(lambda x: 'red' if x.mapq==0 else fill_color ,1)
    df['span'] = df.apply(lambda x: str(x.start)+':'+str(x.end),1)

    if ax==None:
        fig,ax = plt.subplots(1,1,figsize=(15,3))
    from matplotlib.collections import PatchCollection
    patches=[]
    for i,r in df.iterrows():
        rect = plt.Rectangle((r.x, r.y), r.length, h,
                                alpha=.6, linewidth=.5,
                                edgecolor='black', facecolor=r.color)
        patches.append(rect)

    #cmap = ListedColormap(list(df.color))
    ax.add_collection(PatchCollection(patches, match_original=True))
    ax.set_ylim(ystart,yend)
    ax.set_xlim(xstart, xend)
    plt.yticks([])
    plt.tight_layout()
    return

def plot_features(rec, ax, rows=3, xstart=0, xend=30000):

    h=1
    df = tools.records_to_dataframe([rec])
    df = df[(df.feat_type!='region') & (df['feat_type']!='source')]
    df = df[(df.start>xstart) & (df.end<xend)]
    df['length'] = df.end-df.start
    y = list(range(1,rows)) * len(df)
    df['y'] = y[:len(df)]
    df['color'] = 'blue'
    df = df.fillna('')
    #print (df)

    from matplotlib.collections import PatchCollection
    import matplotlib.patches as mpatches

    patches=[]
    for i,r in df.iterrows():
        if r.strand == 1:
            x = r.start
            dx = r.length
        else:
            x = r.end
            dx = -r.length
        arrow = mpatches.Arrow(x, r.y, dx, 0, alpha=.7, width=.3,
                               edgecolor='black')
        txt = ax.text(r.start, r.y-h/2, r.gene, size=16)
        patches.append(arrow)

    ax.add_collection(PatchCollection(patches, match_original=True))
    ax.set_xlim(xstart, xend)
    ax.set_ylim(.4,rows-.5)
    plt.yticks([])
    plt.tight_layout()

    def onclick(event):
        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
              ('double' if event.dblclick else 'single', event.button,
               event.x, event.y, event.xdata, event.ydata))
        ax.text(event.x, event.y, 'HI!')
        ax.figure.canvas.draw()
    #cid = ax.figure.canvas.mpl_connect('button_press_event', onclick)
    return
