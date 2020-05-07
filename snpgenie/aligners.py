"""
    Aligner methods for bacterial genomics.
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

from __future__ import absolute_import, print_function
import sys, os, string, types, re
import platform
import shutil, glob, collections
import itertools
import subprocess
import numpy as np
import pandas as pd
from . import tools

def build_bwa_index(fastafile, path=None):
    """Build a bwa index"""

    bwacmd = tools.get_cmd('bwa')
    cmd = '{b} index {i}'.format(b=bwacmd,i=fastafile)
    subprocess.check_output(cmd, shell=True)
    print (cmd)
    return

def bwa_align(file1, file2, idx, out, threads=4, overwrite=False, filter=None):
    """Align reads to a reference with bwa.
    Args:
        file1, file2: fastq files
        idx: bwa index name
        out: output bam file name
    """

    bwacmd = tools.get_cmd('bwa')
    samtoolscmd = tools.get_cmd('samtools')
    if file2 == None:
        file2=''
    cmd = '{b} mem -M -t {t} {i} {f1} {f2} | {s} view -F 0x04 -bt - | {s} sort -o {o}'.format(b=bwacmd,i=idx,s=samtoolscmd,
                                                                                      f1=file1,f2=file2,o=out,t=threads)
    if not os.path.exists(out) or overwrite == True:
        print (cmd)
        tmp = subprocess.check_output(cmd, shell=True)
    return

def build_bowtie_index(fastafile, path=None):
    """Build a bowtie index
    Args:
        fastafile: file input
        path: folder to place index files
    """

    name = os.path.splitext(os.path.basename(fastafile))[0]
    name = os.path.join(path, name)
    if not os.path.exists(path):
        os.makedirs(path)
    cmd = 'bowtie-build -f %s %s' %(fastafile, name)
    try:
        result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        print (str(e.output))
        return
    print ('built bowtie index for %s' %fastafile)
    return

def build_subread_index(fastafile, path):
    """Build an index for subread"""

    name = os.path.splitext(fastafile)[0]
    cmd = 'subread-buildindex -o %s %s' %(name,fastafile)
    try:
        result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        #print (str(e.output))
        return
    exts = ['.00.b.array','.00.b.tab','.files','.reads']
    files = [name+i for i in exts]
    utils.move_files(files, path)
    return

def bowtie_align(infile, ref, outfile=None, remaining=None, threads=2, verbose=True):
    """Map reads using bowtie"""

    label = os.path.splitext(os.path.basename(infile))[0]
    outpath = os.path.dirname(os.path.abspath(infile))
    if outfile == None:
        outfile = label+'_'+ref+'_bowtie.sam'

    if BOWTIE_INDEXES == None:
        print ('aligners.BOWTIE_INDEXES variable not set')
        return
    os.environ["BOWTIE_INDEXES"] = BOWTIE_INDEXES
    #print (BOWTIE_INDEXES)
    params = BOWTIE_PARAMS
    if remaining == None:
        remaining = os.path.join(outpath, label+'_r.fa')
    cmd = 'bowtie -f -p %s -S %s --un %s %s %s > %s' %(threads,params,remaining,ref,infile,outfile)
    if verbose == True:
        print (cmd)
    try:
        result = subprocess.check_output(cmd, shell=True, executable='/bin/bash',
                                         stderr= subprocess.STDOUT)
        if verbose == True:
            print (result.decode())
    except subprocess.CalledProcessError as e:
        print (str(e.output))
    return remaining
