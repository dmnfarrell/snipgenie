#!/usr/bin/env python

"""
    pathogenie methods for pipeline and cmd line tool.
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
import sys,os,subprocess,glob,re
import platform
import urllib, hashlib, shutil
import tempfile
import pandas as pd
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna
from . import tools, aligners

tempdir = tempfile.gettempdir()
home = os.path.expanduser("~")
config_path = os.path.join(home,'.config/pathogenie')
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'data')
dbdir = os.path.join(config_path, 'db')
customdbdir = os.path.join(config_path, 'custom')
ref_genome = os.path.join(datadir,'Mbovis_AF212297.fa')

def fetch_binaries():
    """Get windows binaries -- windows only"""

    url = "https://github.com/dmnfarrell/btbgenie/pathogenie/raw/master/win_binaries/"
    path = os.path.join(config_path, 'binaries')
    os.makedirs(path, exist_ok=True)
    names = ['blastn.exe','makeblastdb.exe',
            'bcftools.exe','bwa.exe','samtools.exe',
            'msys-2.0.dll','msys-bz2-1.dll','msys-lzma-5.dll','msys-ncursesw6.dll','msys-z.dll']
    for n in names:
        filename = os.path.join(path,n)
        if os.path.exists(filename):
            continue
        link = os.path.join(url,n)
        print (filename,link)
        urllib.request.urlretrieve(link, filename)
    return

def get_sample_names(filenames):
    """Get sample pairs from list of fastq files."""

    res = []
    cols = ['name','sample','filename']
    for filename in filenames:
        name = os.path.basename(filename).split('.')[0]
        sample = name.split('_R')[0]
        x = [name, sample, filename]
        res.append(x)

    df = pd.DataFrame(res, columns=cols)
    df['pair'] = df.groupby('sample').cumcount()+1
    df = df.sort_values(['name','sample','pair']).reset_index(drop=True)
    return df

def align_reads(samples, idx, outdir='mapped', callback=None, **kwargs):
    """
    Align multiple files. Requires a dataframe with a 'sample' column to indicate
    paired files grouping.
    """

    bwacmd = 'bwa'
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    new = []
    for name,df in samples.groupby('sample'):
        print (name)
        if callback != None:
            callback('aligning %s' %name)
        files = list(df.filename)
        #print (files)
        out = os.path.join(outdir,name+'.bam')
        print (out)
        aligners.bwa_align(files[0],files[1], idx=idx, out=out, **kwargs)
        cmd = 'samtools index {o}'.format(o=out)
        subprocess.check_output(cmd,shell=True)
        index = df.index
        samples.loc[index,'bam_file'] = out
        if callback != None:
            callback(out)
    return samples

def variant_calling(bam_files, ref, outpath, callback=None, **kwargs):
    """Call variants with bcftools"""

    bcftoolscmd = 'bcftools'
    if platform.system() == 'Windows':
        bcftoolscmd = tools.win_binary_path('bcftools.exe')

    rawbcf = os.path.join(outpath,'raw.bcf')
    cmd = '{c} mpileup -O b -o {o} -f {r} {b}'.format(r=ref, b=bam_files, o=rawbcf,c=bcftoolscmd)
    print (cmd)
    if callback != None:
        callback(cmd)
    subprocess.check_output(cmd,shell=True)
    #find snps
    vcfout = os.path.join(outpath,'calls.vcf')
    cmd = '{c} call --ploidy 1 -m -v -o {v} {raw}'.format(v=vcfout,raw=rawbcf,c=bcftoolscmd)
    if callback != None:
        callback(cmd)
    print (cmd)
    subprocess.check_output(cmd,shell=True)
    #filter the calls
    vcfout = os.path.join(outpath,'calls.vcf')
    final = os.path.join(outpath,'filtered')
    cmd = 'vcftools --vcf {i} --minQ 20 --recode --recode-INFO-all --out {o}'.format(i=vcfout,o=final,c=bcftoolscmd)
    print (cmd)
    tmp = subprocess.check_output(cmd,shell=True)
    if callback != None:
        callback(tmp)
    return final
