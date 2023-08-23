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
import matplotlib as mpl
import matplotlib.colors as colors
from . import tools

def make_legend(fig, colormap, loc=(1.05, .6), title='',fontsize=12):
    """Make a figure legend wth provided color mapping"""

    import matplotlib.patches as mpatches
    pts=[]
    for c in colormap:
        pts.append(mpatches.Patch(color=colormap[c],label=c))
    fig.legend(handles=pts,bbox_to_anchor=loc,fontsize=fontsize,title=title)
    return pts

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

def random_colors(n=10, seed=1):
    """Generate random hex colors as list of length n."""

    import random
    random.seed(seed)
    clrs=[]
    for i in range(n):
        r = lambda: random.randint(0,255)
        c='#%02X%02X%02X' % (r(),r(),r())
        clrs.append(c)
    return clrs

def random_cmap(n, seed=10):

    rcolors = random_colors(n,seed)
    c=range(n)
    return dict(zip(c, rcolors))

def gen_colors(cmap,n,reverse=False):
    '''Generates n distinct color from a given colormap.
    Args:
        cmap(str): The name of the colormap you want to use.
            Refer https://matplotlib.org/stable/tutorials/colors/colormaps.html to choose
            Suggestions:
            For Metallicity in Astrophysics: Use coolwarm, bwr, seismic in reverse
            For distinct objects: Use gnuplot, brg, jet,turbo.
        n(int): Number of colors you want from the cmap you entered.
        reverse(bool): False by default. Set it to True if you want the cmap result to be reversed.
    Returns:
        colorlist(list): A list with hex values of colors.

    Taken from the mycolorpy package by binodbhttr
    see also https://matplotlib.org/stable/tutorials/colors/colormaps.html
    '''

    c_map = plt.cm.get_cmap(str(cmap)) # select the desired cmap
    arr=np.linspace(0,1,n) #create a list with numbers from 0 to 1 with n items
    colorlist=list()
    for c in arr:
        rgba=c_map(c) #select the rgba value of the cmap at point c which is a number between 0 to 1
        clr=colors.rgb2hex(rgba) #convert to hex
        colorlist.append(str(clr)) # create a list of these colors

    if reverse==True:
        colorlist.reverse()
    return colorlist

def show_colors(colors):
    """display a list of colors"""

    plt.figure(figsize=(6,1))
    plt.bar(range(len(colors)),height=1,color=colors,width=1)
    plt.axis('off')
    return

def get_color_mapping(df, col, cmap=None, seed=1):
    """Get random color map for categorcical dataframe column"""

    c = df[col].unique()
    if cmap == None:
        rcolors = random_colors(len(c),seed)
    else:
        cmap = mpl.cm.get_cmap(cmap)
        rcolors = [cmap(i) for i in range(len(c))]
    colormap = dict(zip(c, rcolors))
    newcolors =  [colormap[i] if i in colormap else 'Black' for i in df[col]]
    return newcolors, colormap

def draw_pie(vals, xpos, ypos, colors, size=500, ax=None):
    """Draw a pie at a specific position on an mpl axis.
    Used to draw spatial pie charts on maps.
    Args:
        vals: values for pie
        xpos: x coord
        ypos: y coord
        colors: colors of values
        size: size of pie chart
    """

    cumsum = np.cumsum(vals)
    cumsum = cumsum / cumsum[-1]
    pie = [0] + cumsum.tolist()

    #colors = ["blue", "red", "yellow"]
    for i, (r1, r2) in enumerate(zip(pie[:-1], pie[1:])):
        angles = np.linspace(2 * np.pi * r1, 2 * np.pi * r2)
        x = [0] + np.cos(angles).tolist()
        y = [0] + np.sin(angles).tolist()
        xy = np.column_stack([x, y])
        ax.scatter([xpos], [ypos], marker=xy, s=size, color=colors[i], alpha=.9)

    return ax

def create_grid(gdf=None, bounds=None, n_cells=10, overlap=False, crs="EPSG:29902"):
    """Create square grid that covers a geodataframe area
    or a fixed boundary with x-y coords
    returns: a GeoDataFrame of grid polygons
    """

    import geopandas as gpd
    import shapely

    if bounds != None:
        xmin, ymin, xmax, ymax= bounds
    else:
        xmin, ymin, xmax, ymax= gdf.total_bounds

    # get cell size
    cell_size = (xmax-xmin)/n_cells
    # create the cells in a loop
    grid_cells = []
    for x0 in np.arange(xmin, xmax+cell_size, cell_size ):
        for y0 in np.arange(ymin, ymax+cell_size, cell_size):
            x1 = x0-cell_size
            y1 = y0+cell_size
            poly = shapely.geometry.box(x0, y0, x1, y1)
            #print (gdf.overlay(poly, how='intersection'))
            grid_cells.append( poly )

    cells = gpd.GeoDataFrame(grid_cells, columns=['geometry'],
                                     crs=crs)
    if overlap == True:
        cols = ['grid_id','geometry','grid_area']
        cells = cells.sjoin(gdf, how='inner').drop_duplicates('geometry')
    return cells

def create_hex_grid(gdf=None, bounds=None, n_cells=10, overlap=False, crs="EPSG:29902"):
    """Hexagonal grid over geometry.
    See https://sabrinadchan.github.io/data-blog/building-a-hexagonal-cartogram.html
    """

    from shapely.geometry import Polygon
    import geopandas as gpd
    if bounds != None:
        xmin, ymin, xmax, ymax= bounds
    else:
        xmin, ymin, xmax, ymax= gdf.total_bounds

    unit = (xmax-xmin)/n_cells
    a = np.sin(np.pi / 3)
    cols = np.arange(np.floor(xmin), np.ceil(xmax), 3 * unit)
    rows = np.arange(np.floor(ymin) / a, np.ceil(ymax) / a, unit)

    #print (len(cols))
    hexagons = []
    for x in cols:
      for i, y in enumerate(rows):
        if (i % 2 == 0):
          x0 = x
        else:
          x0 = x + 1.5 * unit

        hexagons.append(Polygon([
          (x0, y * a),
          (x0 + unit, y * a),
          (x0 + (1.5 * unit), (y + unit) * a),
          (x0 + unit, (y + (2 * unit)) * a),
          (x0, (y + (2 * unit)) * a),
          (x0 - (0.5 * unit), (y + unit) * a),
        ]))

    grid = gpd.GeoDataFrame({'geometry': hexagons},crs=crs)
    grid["grid_area"] = grid.area
    grid = grid.reset_index().rename(columns={"index": "grid_id"})
    if overlap == True:
        cols = ['grid_id','geometry','grid_area']
        grid = grid.sjoin(gdf, how='inner').drop_duplicates('geometry')[cols]
    return grid


# below are sequence plotting related functions

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

def get_chrom_from_bam(bam_file):
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

def display_igv(url='http://localhost:8888/files/', ref_fasta='', bams=[],
                gff_file=None, vcf_file=None):
    """
    Display IGV tracks in jupyter, requires the igv_jupyterlab package.
    Example usage:
        bams = {'24':'results/mapped/24-MBovis.bam'}
        igv=display_igv(url='http://localhost:8888/files/', ref_fasta='Mbovis_AF212297.fa',
            gff_file='results/Mbovis_AF212297.gb.gff',
            vcf_file='results/filtered.vcf.gz', bams=bams)
    """

    from igv_jupyterlab import IGV

    track_list = []
    if gff_file != None:
        track_list += [{"name": "Mbovis",
                        "url": url+gff_file,
                        "format": "gff",
                        "type": "annotation",
                        "height":120,
                        "indexed": False }
                     ]

    colors = random_colors(len(bams))
    i=0
    for b in bams:
        d = {"name": b,
            "url":url+bams[b],
            "type": "alignment",
             "displayMode":"SQUISHED",
             "height":110,
             "removable":True,
             "color":colors[i],
            "indexed": True }
        track_list.append(d)
        i+=1

    if vcf_file != None:
        v = IGV.create_track(
            name="VCF",
            url=url+vcf_file,
            #format="vcf",
            type="variant",
            indexed=True
        )
        track_list.append(v)

    genome = IGV.create_genome(
        name="Mbovis",
        fasta_url=url+ref_fasta,
        index_url=url+ref_fasta+'.fai',
        tracks=track_list
    )

    igv = IGV(genome=genome)
    #igv.locus="LT708304.1:2283938-2285249"

    return igv

bams={'cat':'wicklow_results/mapped/cat-003488.bam',
      'cat2018':'wicklow_results/mapped/TB18-001924.bam',
      'cat2020':'wicklow_results/mapped/TB20-003486.bam',
      '24':'wicklow_results/mapped/24-MBovis.bam'
      }
