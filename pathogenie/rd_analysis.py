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
import platform
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Phylo, AlignIO
import numpy as np
import pandas as pd
from gzip import open as gzopen
import tempfile
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
    df=RD.set_index('RD_name')
    if names!= None:
        df=df.loc[names]
    seqs=[]
    for name, row in df.iterrows():
        #print (name,row.Start, row.Stop, row.Rv)
        from pyfaidx import Fasta
        rg = Fasta('../MTB-H37Rv.fna')
        sseq = rg['NC_000962.3'][row.Start:row.Stop].seq
        #refname = '%s.fa' %name
        seqs.append(SeqRecord(Seq(sseq),id=name))
    SeqIO.write(seqs, 'RD.fa', 'fasta')
    aligners.build_bwa_index('RD.fa')

def align_regions(df, path):
    """Align reads to regions of difference"""

    from io import StringIO
    from pyfaidx import Fasta
    ref = 'RD.fa'
    rg = Fasta('../MTB-H37Rv.fna')
    res = []
    for i,g in df.groupby('sample'):
        out=os.path.join(path,i+'.bam')
        f1 = g.iloc[0].filename; f2 = g.iloc[1].filename
        if not os.path.exists(out):
            aligners.bwa_align(f1, f2, ref, out, threads=4, overwrite=False)

        cmd = 'zcat %s | paste - - - - | cut -f2 | wc -c' %f1
        tmp = subprocess.check_output(cmd,shell=True)
        avdepth = int(tmp)*2/len(rg)
        print (avdepth)
        cmd = 'samtools coverage --min-BQ 1 %s' %out
        tmp = subprocess.check_output(cmd,shell=True)
        s = pd.read_csv(StringIO(tmp.decode()),sep='\t')
        s['name'] = i
        #print (s)
        s['ratio'] = s.meandepth/avdepth
        res.append(s)
    res = pd.concat(res)
    return res
