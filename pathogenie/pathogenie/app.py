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
sequencedir = os.path.join(config_path, 'genome')
ref_genome = os.path.join(datadir,'Mbovis_AF212297.fa')
ref_gff = os.path.join(datadir,'Mbovis_AF212297.gff')

def copy_ref_genomes():
    """Copy default ref genome files to config dir"""

    path = sequencedir
    if not os.path.exists(path):
        os.makedirs(path,exist_ok=True)
    dest = os.path.join(path, os.path.basename(ref_genome))
    shutil.copy(ref_genome, dest)
    return

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

def get_sample_names(filenames, sep='-'):
    """Get sample pairs from list of fastq files."""

    res = []
    cols = ['name','sample','filename']
    for filename in filenames:
        name = os.path.basename(filename).split('.')[0]
        sample = name.split(sep)[0]
        #if we can't get sample name try another delimeter?
        if name == sample:
            sample = name.split('_')[0]
        x = [name, sample, filename]
        res.append(x)

    df = pd.DataFrame(res, columns=cols)
    df = df.sort_values(['name','sample']).reset_index(drop=True)
    df['pair'] = df.groupby('sample').cumcount()+1
    #df = df.sort_values(['name','sample','pair']).reset_index(drop=True)
    return df

def align_reads(samples, idx, outdir='mapped', callback=None, **kwargs):
    """
    Align multiple files. Requires a dataframe with a 'sample' column to indicate
    paired files grouping. If a trimmed column is present these files will align_reads
    instead of the raw ones.
    """

    bwacmd = 'bwa'
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    new = []
    for name,df in samples.groupby('sample'):
        print (name)
        if callback != None:
            callback('aligning %s' %name)
        if 'trimmed' in df.columns:
            files = list(df.trimmed)
            callback('using trimmed')
        else:
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

def worker(region,out,bam_files):
    """Run bcftools for single region."""

    cmd = 'bcftools mpileup -r {reg} -O b -o {o} -f {r} {b}'.format(r=ref, reg=region, b=bam_files, o=out)
    #print (cmd)
    subprocess.check_output(cmd, shell=True)
    cmd = 'bcftools index {o}'.format(o=out)
    subprocess.check_output(cmd, shell=True)
    return

def mpileup_parallel(bam_files, ref, outpath, threads=4, callback=None):
    """Run mpileup in parallel over multiple regions, then concat vcf files.
    Assumes alignment to a bacterial reference with a single chromosome."""

    import multiprocessing as mp
    bam_files = ' '.join(bam_files)
    ref = app.ref_genome
    rawbcf = os.path.join(outpath,'raw.bcf')
    tmpdir = 'tmp'
    chr = 'NC_002945.4'
    length = tools.get_fasta_length(ref)
    #print (length)
    #find regions
    bsize = int(length/(threads-1))
    x = np.linspace(1,length,threads,dtype=int)
    blocks=[]
    for i in range(len(x)):
        if i < len(x)-1:
            blocks.append((x[i],x[i+1]-1))
    #print (blocks, bsize)

    pool = mp.Pool(threads)
    res = []
    outfiles = []
    st = time.time()

    for start,end in blocks:
        #end = start+bsize
        print (start, end)
        region = '{c}:{s}-{e}'.format(c=chr,s=start,e=end)
        out = '{o}/{s}.bcf'.format(o=tmpdir,s=start)
        f = pool.apply_async(worker, [region,out,bam_files])
        res.append(f)
        outfiles.append(out)

    pool.close()
    pool.join()
    t=time.time()-st
    print ('took %s seconds' %str(round(t,3)))

    #concat files
    cmd = 'bcftools concat {i} -O b -o {o}'.format(i=' '.join(outfiles),o=rawbcf)
    print (cmd)
    subprocess.check_output(cmd, shell=True)
    #remove temp files
    for f in outfiles:
        os.remove(f)
    return rawbcf

def variant_calling(bam_files, ref, outpath, sample_file=None, callback=None, overwrite=False, **kwargs):
    """Call variants with bcftools"""

    bam_files = ' '.join(bam_files)
    rawbcf = os.path.join(outpath,'raw.bcf')
    cmd = 'bcftools mpileup -O b -o {o} -f {r} {b}'.format(r=ref, b=bam_files, o=rawbcf)

    if callback != None:
        callback(cmd)
    if not os.path.exists(rawbcf) or overwrite == True:
        print (cmd)
        subprocess.check_output(cmd,shell=True)
    #find snps
    vcfout = os.path.join(outpath,'calls.vcf')
    cmd = 'bcftools call --ploidy 1 -m -v -o {v} {raw}'.format(v=vcfout,raw=rawbcf)
    if callback != None:
        callback(cmd)
    print (cmd)
    subprocess.check_output(cmd,shell=True)
    #rename samples
    if sample_file != None:
        cmd = 'bcftools reheader --samples {s} -o {v} {v}'.format(v=vcfout,s=sample_file)
        print(cmd)
        tmp = subprocess.check_output(cmd,shell=True)
    final = os.path.join(outpath,'filtered')
    #cmd = 'vcftools --vcf {i} --minQ 20 --recode --recode-INFO-all --out {o}'.format(i=vcfout,o=final)
    cmd = 'bcftools filter -e "QUAL<40" -o {o}.vcf.gz -O z {i}'.format(i=vcfout,o=final)
    print (cmd)
    tmp = subprocess.check_output(cmd,shell=True)
    if callback != None:
        callback(tmp)
    return vcfout

def create_bam_labels(filenames):

    names = [os.path.basename(i).split('.')[0] for i in filenames]
    with open('samples.txt','w+') as file:
        for s in zip(bam_files,names):
            file.write('%s %s\n' %(s[0],s[1]))
    return

def fasta_alignment_from_vcf(vcf_file, ref, callback=None):
    """Get a fasta alignment for all snp sites in a multi sample
    vcf file, including the reference sequence"""

    from pyfaidx import Fasta
    from pyfaidx import FastaVariant
    #index vcf
    cmd = 'tabix -p vcf -f {i}'.format(i=vcf_file)
    tmp = subprocess.check_output(cmd,shell=True)
    #get samples?
    import vcf
    vcf_reader = vcf.Reader(open(vcf_file, 'rb'))
    samples = vcf_reader.samples
    print ('%s samples' %len(samples))
    result = []

    #reference sequence
    reference = Fasta(ref)
    chrom = list(reference.keys())[0]

    #get the set of all sites first
    sites=[]
    for sample in samples:
        #print (sample)
        variant = FastaVariant(ref, vcf_file,
                                 sample=sample, het=True, hom=True)
        pos = list(variant[chrom].variant_sites)
        sites.extend(pos)
        #print (sample)
        #print (pos[:20])
    sites = sorted(set(sites))
    n=550
    #sites=sites[:n]
    print ('using %s sites' %len(sites))
    if callback != None:
        callback('using %s sites' %len(sites))
    #get reference sequence for site positions
    refseq=[]
    for p in sites:
        refseq.append(reference[chrom][p-1].seq)
    refseq = ''.join(refseq)
    #print (refseq)
    refrec = SeqRecord(Seq(refseq),id='ref')
    result.append(refrec)

    #iterate over variants in each sample
    for sample in samples:
        seq=[]
        variant = FastaVariant(ref, vcf_file,
                                 sample=sample, het=True, hom=True)
        #for p in variant[chrom].variant_sites:
        for p in sites:
            rec = variant[chrom][p-1:p]
            #print (p,rec)
            seq.append(rec.seq)
        seq = ''.join(seq)
        #print (seq)
        seqrec = SeqRecord(Seq(seq),id=sample)
        result.append(seqrec)
    return result
