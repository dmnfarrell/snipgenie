"""
    Region of difference analysis for MTBC isolates.
    see https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3213-1
    Created March 2020
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

import sys,os,io,subprocess,glob,shutil,re,random
import numpy as np
import pandas as pd
import pylab as plt
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from . import tools, aligners

tempdir = tempfile.gettempdir()
home = os.path.expanduser("~")
config_path = os.path.join(home,'.config','pathogenie')
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'data')
mtbref = os.path.join(datadir, 'MTB-H37Rv.fa')

def create_rd_index(names=None):
    """Get RD region sequence from reference and make bwa index"""

    RD = pd.read_csv(os.path.join(datadir,'RD.csv'))
    df = RD.set_index('RD_name')
    if names!= None:
        df=df.loc[names]
    seqs=[]
    for name, row in df.iterrows():
        #print (name,row.Start, row.Stop, row.Rv)
        from pyfaidx import Fasta
        rg = Fasta(mtbref)
        sseq = rg['NC_000962.3'][row.Start:row.Stop].seq
        #refname = '%s.fa' %name
        seqs.append(SeqRecord(Seq(sseq),id=name))
    SeqIO.write(seqs, 'RD.fa', 'fasta')
    aligners.build_bwa_index('RD.fa')

def run_samples(df, path, threads=4):
    """Run a set of samples
    Args:
        df: a samples dataframe from snpgenie
        path: folder with raw reads
    """

    res = []
    df=df.fillna('')
    for i,r in df.iterrows():
        name = r['sample']
        print (name)
        f1 = r.filename1
        f2 = r.filename2
        s = find_regions(f1, f2, path, name, threads)
        res.append(s)
    res = pd.concat(res)
    return res

def get_average_depth(ref, f1):
    from pyfaidx import Fasta
    rg = Fasta(ref)
    k = list(rg.keys())[0]
    cmd = 'zcat "%s" | paste - - - - | cut -f2 | wc -c' %f1
    tmp = subprocess.check_output(cmd,shell=True)
    avdepth = int(tmp)*2/len(rg[k])
    return avdepth

def find_regions(f1, f2, path, name, threads=4, avdepth=None):
    """Align reads to regions of difference and get coverage stats.
    Args:
        f1: first filename
        path: folder with raw reads
    Returns:
        list of mapping values with depths
    """

    from io import StringIO
    from pyfaidx import Fasta
    ref = 'RD.fa'
    #rg = Fasta(mtbref)
    #k = list(rg.keys())[0]
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
    out = os.path.join(path,name+'.bam')
    if not os.path.exists(out):
        aligners.bwa_align(f1, f2, ref, out, threads=threads, overwrite=False)
    #get the average sequencing depth
    if avdepth == None:
        avdepth = get_average_depth(mtbref, f1)
    #print (avdepth)
    cmd = 'samtools coverage --min-BQ 0 %s' %out
    tmp = subprocess.check_output(cmd,shell=True)
    s = pd.read_csv(StringIO(tmp.decode()),sep='\t')
    s['name'] = name
    s['ratio'] = s.meandepth/avdepth
    return s

def get_matrix(res, cutoff=0.15):
    """Get presence/absence matrix for RDs. Cutoff is ratio of reads to region vs
    mean depth. and is arbitrary."""

    X = pd.pivot_table(res,index='name',columns=['#rname'],values='ratio')
    X=X.clip(lower=cutoff).replace(cutoff,0)
    X=X.clip(upper=cutoff).replace(cutoff,1)
    #invert so 1 means RD is present
    X = X.replace({0:1, 1:0})
    X=X.sort_values(by=X.columns[0])
    return X

def apply_rules(x):
    """Identify isolate using RD rules"""

    if x.RD239 == 0:
        return 'L1'
    elif x.RD105 == 0:
        return 'L2'
    elif x.RD4 == 1:
        if (x.RD1mic == 0):
            return 'Microti'
        elif (x.RD12bov == 0 or x.RD1bcg == 0 or x.RD2bcg == 0):
            return 'Caprae'
    elif x.RD4 == 0:
        if x.RD1bcg==1 and x.RD2bcg==1 and x.RD12bov == 0:
            return 'Bovis'
        elif x.RD1bcg==0 and x.RD2bcg==1:
            return 'BCG (Moreau)'
        elif x.RD1bcg==0 and x.RD2bcg==0:
            return 'BCG (Merieux)'
    elif x.RD711 == 0:
        return 'Africanum'
    else:
        return 'Unknown'

def get_coverage(bam_file, chr, start, end, ref):
    """Get coverage of a region. Returns a dataframe"""

    cmd = 'samtools mpileup {b} --min-MQ 10 -f {r} -r {c}:{s}-{e}'.format(c=chr,s=start,e=end,b=bam_file,r=ref)
    #print(cmd)
    temp = subprocess.check_output(cmd, shell=True)
    df=pd.read_csv(io.BytesIO(temp), sep='\t', names=['chr','pos','base','coverage','q','c'])
    if len(df)==0:
        df['pos'] = range(start,end)
        df['coverage']=0
    return df

def show_rd_coverage(samples, chr, start, end, ref_fasta, ref_gb, margin=None,
                    labelcol=None,title=None, colors={}):
    """
    Plot read coverage in specific regions for multiple samples. useful for
    RD region comparison.
    """

    from dna_features_viewer import GraphicFeature, GraphicRecord
    from dna_features_viewer import BiopythonTranslator
    from matplotlib.gridspec import GridSpec
    import matplotlib.ticker as ticker

    rec = list(SeqIO.parse(ref_gb,format='gb'))[0]
    rec.features = [f for f in rec.features if f.type!='gene']
    rd = rec[start:end]
    graphic_record = BiopythonTranslator().translate_record(rd)

    fig = plt.figure(figsize=(22,8))
    gs = GridSpec(len(samples)+3, 1, figure=fig)
    ax1=fig.add_subplot(gs[:2,0])
    graphic_record.plot(ax=ax1,with_ruler=False)
    i=3
    if margin==None:
        margin=(end-start)*.25

    for n,r in samples.iterrows():
        name=r['sample']
        ax=fig.add_subplot(gs[i,0])
        df = get_coverage(r.bam_file,chr,start-margin,end+margin,ref_fasta)
        df=df.set_index('pos')
        #print (df)
        #bins=range(0,max(df.coverage),int(max(df.coverage)/5))
        #df['binned']=np.searchsorted(bins, df.coverage.values)
        if name in colors:
            color=colors[name]
        else:
            color='gray'
        df['coverage']=df.coverage.replace(0,np.nan)
        df.plot(y='coverage',ax=ax,kind='area',color=color,legend=False)
        #print (df)
        if labelcol != None:
            label = r[labelcol]
        else:
            label = name
        ax.text(-.12,.5,label,color='blue',transform=ax.transAxes,fontsize=12)
        ax.set_yticklabels([])
        ax.set_xlim(start-margin, end+margin)
        ax.axvline(x=start);ax.axvline(x=end)
        i+=1

    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    plt.subplots_adjust(left=.3,right=.9,wspace=0, hspace=0)
    if title==None:
        title='%s:%s-%s' %(chr,start,end)
    fig.suptitle(title,fontsize=25)
    return
