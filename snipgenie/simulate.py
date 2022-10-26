#!/usr/bin/env python

"""
    Simulate reads
    Created Sep 2022
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

import sys, os, io
import numpy as np
import string
import subprocess
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO, AlignIO
from .qt import *
from . import tools, widgets

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
config_path = os.path.join(home,'.config','snipgenie')
bin_path = os.path.join(config_path, 'binaries')

def run_phastsim(path, ref, newick):
    """Run phastsim """

    refseq = SeqIO.read(ref, 'fasta')
    scale = (10/len(refseq))
    outpath = os.path.join(path, 'phastsim_output')
    cmd = 'phastSim --outpath {o} --outputFile 1 --seed 1 --createFasta' \
             ' --createPhylip --treeFile {n}' \
             ' --scale {s} --invariable .1 --alpha 1.0 --omegaAlpha 1.0' \
             ' --reference {r}'.format(o=outpath,s=scale,r=ref,n=newick)
    print (cmd)
    subprocess.check_output(cmd, shell=True)
    return

def artificial_fastq_generator(ref, outfile, cmp=100):
    """Generate reads from reference"""

    jarfile = os.path.join(bin_path,'ArtificialFastqGenerator.jar')
    cmd = 'java -jar {j} -O {o} -R {r} -S ">temp" -RL 150 -CMP {cmp}'\
        ' -CSD 0.2 -SE true'.format(j=jarfile, r=ref, o=outfile,cmp=cmp)
         #-URQS true -F1 {f1} -F2 {f2}
    print (cmd)
    subprocess.check_output(cmd, shell=True)
    return

def generate_fastqs(infile, outpath, reads=1e5, overwrite=False):
    """Make multiple fastqs"""

    from joblib import Parallel, delayed
    import multiprocessing, time
    num_cores = 4

    simrecs = list(SeqIO.parse(infile,'fasta'))
    def my_func(rec):
        from tempfile import mkstemp
        x,tmp = mkstemp()
        SeqIO.write(SeqRecord(rec.seq,id='temp'), tmp, 'fasta')
        out = os.path.join(outpath,rec.id)
        if os.path.exists(out+'.1.fastq.gz') and overwrite == False:
            print ('found %s' %out)
            return
        cmp = reads*2*150/len(rec.seq)
        print (cmp)
        artificial_fastq_generator(tmp, out, cmp=cmp)
        cmd = 'pigz %s/*.fastq' %outpath
        subprocess.check_output(cmd, shell=True)

    st = time.time()
    #Parallel(n_jobs=num_cores)(delayed(my_func)(i) for i in simrecs)
    #print (time.time()-st)
    for i in simrecs:
        my_func(i)
