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

from __future__ import print_function
import sys,os,subprocess,glob,shutil,re,random
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

def find_regions(df, path, threads=4, callback=None):
    """Align reads to regions of difference and get coverage stats."""

    from io import StringIO
    from pyfaidx import Fasta
    ref = 'RD.fa'
    rg = Fasta(mtbref)
    k = list(rg.keys())[0]
    res = []
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
    for i,g in df.groupby('sample'):
        out = os.path.join(path,i+'.bam')
        f1 = g.iloc[0].filename; f2 = g.iloc[1].filename
        if not os.path.exists(out):
            aligners.bwa_align(f1, f2, ref, out, threads=threads, overwrite=False)
        #get the average sequencing depth
        cmd = 'zcat %s | paste - - - - | cut -f2 | wc -c' %f1
        tmp = subprocess.check_output(cmd,shell=True)
        print (rg)
        avdepth = int(tmp)*2/len(rg[k])
        print (avdepth)
        cmd = 'samtools coverage --min-BQ 0 %s' %out
        tmp = subprocess.check_output(cmd,shell=True)
        s = pd.read_csv(StringIO(tmp.decode()),sep='\t')
        s['name'] = i
        #print (s)
        s['ratio'] = s.meandepth/avdepth
        res.append(s)
        if callback != None:
            callback(i)
    res = pd.concat(res)
    return res

def get_matrix(res, cutoff=0.15):
    """Get presence/absence matrix for RDs"""

    X = pd.pivot_table(res,index='name',columns=['#rname'],values='ratio')
    X=X.clip(lower=cutoff).replace(cutoff,0)
    X=X.clip(upper=cutoff).replace(cutoff,1)
    X=X.sort_values(by=X.columns[0])
    #print (X[:4])
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
