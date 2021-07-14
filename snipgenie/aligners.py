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

import sys, os, string, types, re
import platform, tempfile
import shutil, glob, collections
import itertools
import subprocess
import numpy as np
import pandas as pd
from . import tools

tempdir = tempfile.gettempdir()
home = os.path.expanduser("~")
config_path = os.path.join(home,'.config','snipgenie')
module_path = os.path.dirname(os.path.abspath(__file__))
BOWTIE_INDEXES = os.path.join(config_path, 'genome')
SUBREAD_INDEXES = os.path.join(config_path, 'genome')

def build_bwa_index(fastafile, path=None):
    """Build a bwa index"""

    bwacmd = tools.get_cmd('bwa')
    cmd = '{b} index {i}'.format(b=bwacmd,i=fastafile)
    subprocess.check_output(cmd, shell=True)
    print (cmd)
    return

def bwa_align(file1, file2, idx, out, threads=4, overwrite=False,
              options='', filter=None, unmapped=None):
    """Align reads to a reference with bwa.
    Args:
        file1, file2: fastq files
        idx: bwa index name
        out: output bam file name
        options: extra command line options e.g. -k INT for seed length
        unmapped: path to file for unmapped reads if required
    """

    bwacmd = tools.get_cmd('bwa')
    samtoolscmd = tools.get_cmd('samtools')
    if file2 == None:
        file2=''
    cmd = '{b} mem -M -t {t} {p} {i} "{f1}" "{f2}" | {s} view -bt - | {s} sort -o {o}'.format(
                b=bwacmd,i=idx,s=samtoolscmd,
                f1=file1,f2=file2,o=out,t=threads,p=options)
    if not os.path.exists(out) or overwrite == True:
        print (cmd)
        tmp = subprocess.check_output(cmd, shell=True)

        #write out unmapped reads
        if unmapped != None:
            f1 = os.path.join(unmapped,os.path.basename(file1))
            f2 = os.path.join(unmapped,os.path.basename(file2))
            #cmd = '{s} view -b -f 4 {o} > {u}'.format(s=samtoolscmd,o=out,u=unmapped)
            #cmd = '{s} view -b -f 4 {o} | {s} fastq -1 {f1} -2 {f2} \
            #         -0 /dev/null -s /dev/null -n'.format(s=samtoolscmd,o=out,u=unmapped,f1=f1,f2=f2)

            uf =  os.path.join(unmapped,os.path.basename(out))+'.fasta'
            cmd = '{s} view -b -f 4 {o} | {s} fasta -0 /dev/null {o} > {f}'.format(
                    s=samtoolscmd,o=out,u=unmapped,f=uf)
            #print (cmd)
            tmp = subprocess.check_output(cmd, shell=True)
    return

def build_bowtie_index(fastafile, path=None):
    """Build a bowtie index
    Args:
        fastafile: file input
        path: folder to place index files
    """

    if path == None:
        path = BOWTIE_INDEXES

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
    return name

def bowtie_align(file1, file2, idx, out, remaining=None, threads=2,
                overwrite=False, verbose=True, options='-v 1 --best'):
    """Map reads using bowtie"""

    bowtiecmd = tools.get_cmd('bowtie')
    samtoolscmd = tools.get_cmd('samtools')

    if BOWTIE_INDEXES == None:
        print ('aligners.BOWTIE_INDEXES variable not set')
        return
    os.environ["BOWTIE_INDEXES"] = BOWTIE_INDEXES
    if file2 != None:
        filestr = '-1 {f1} -2 {f2}'.format(f1=file1,f2=file2)
    else:
        filestr = file1
    cmd = '{c} -q -p {t} -S {p} {r} {f} | {s} view -F 0x04 -bt - | {s} sort -o {o}'\
            .format(c=bowtiecmd,t=threads,f=filestr,p=options,r=idx,o=out,s=samtoolscmd)

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

def build_subread_index(fastafile):
    """Build an index for subread"""

    path = SUBREAD_INDEXES
    name = os.path.splitext(fastafile)[0]
    subreadalign = tools.get_cmd('subread-buildindex')
    cmd = '{sc} -o {n} {f}'.format(sc=subreadalign,n=name,f=fastafile)
    print (cmd)
    try:
        result = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print (str(e.output))
        return
    exts = ['.00.b.array','.00.b.tab','.files','.reads']
    files = [name+i for i in exts]
    tools.move_files(files, path)
    return

def subread_align(file1, file2, idx, out, threads=2,
                overwrite=False, verbose=True):

    os.environ["SUBREAD_INDEXES"] = SUBREAD_INDEXES
    idx = os.path.join(SUBREAD_INDEXES, idx)
    samtoolscmd = tools.get_cmd('samtools')
    subreadcmd = tools.get_cmd('subread-align')
    params = '-t 1 --SAMoutput -m 3 -M 2'
    cmd = '{sc} {p} -T {t} -i {i} -r "{f1}" -R "{f2}" | {s} view -F 0x04 -bt - | {s} sort -o {o}'.format(
            sc=subreadcmd,p=params,t=threads,i=idx,f1=file1,f2=file2,s=samtoolscmd,o=out)
    if not os.path.exists(out) or overwrite == True:
        print (cmd)
        result = subprocess.check_output(cmd, shell=True, stderr= subprocess.STDOUT)
    return

def minimap2_align(file, ref, out, threads=4, overwrite=False):
    """Align ONT reads with minimap2"""

    samtoolscmd = tools.get_cmd('samtools')
    minimapcmd = tools.get_cmd('minimap2')
    cmd = '{m} -t {t} -ax map-ont {r} {q} | {s} view -F 0x04 -bt - | {s} sort -o {o}'\
    .format(r=ref,q=query,s=samtoolscmd,m=minimapcmd,o=out,t=threads)
    if not os.path.exists(out) or overwrite == True:
        result = subprocess.check_output(cmd, shell=True, stderr= subprocess.STDOUT)
    return
