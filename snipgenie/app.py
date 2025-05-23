#!/usr/bin/env python

"""
    snipgenie methods for cmd line tool.
    Created Nov 2019
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 3
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warroanty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

import sys,os,subprocess,glob,re
import time, datetime
import platform
import urllib, hashlib, shutil
import tempfile
import pandas as pd
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO, AlignIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
#from Bio.Alphabet import generic_dna
from . import tools, aligners, trees
import multiprocessing as mp

tempdir = tempfile.gettempdir()
home = os.path.expanduser("~")
config_path = os.path.join(home,'.config','snipgenie')
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'data')
sequence_path = os.path.join(config_path, 'genome')
annotation_path = os.path.join(config_path, 'annotation')
mbovis_genome = os.path.join(sequence_path, 'Mbovis_AF212297.fa')
mtb_genome = os.path.join(sequence_path, 'MTB-H37Rv.fa')
mbovis_gb = os.path.join(datadir, 'Mbovis_AF212297.gb')
mtb_gb = os.path.join(datadir, 'MTB-H37Rv.gb')
map_genome = os.path.join(sequence_path, 'MAP-K10.fa')
map_gb = os.path.join(datadir, 'MAP-K10.gb')
msmeg_genome = os.path.join(sequence_path, 'Msmeg-MC2.fa')
msmeg_gb = os.path.join(datadir, 'Msmeg-MC2.gb')
mbovis_mask =  os.path.join(datadir, 'Mbovis_AF212297_mask.bed')
mtb_mask =  os.path.join(datadir, 'MTB-H37Rv_mask.bed')
map_mask =  os.path.join(datadir, 'MAP-K10_mask.bed')
mycoplasmabovis_genome = os.path.join(sequence_path, 'Mpbovis-PG45.fa')
mycoplasmabovis_gb = os.path.join(datadir, 'Mpbovis-PG45.gb')
sarscov2_genome = os.path.join(sequence_path, 'Sars-Cov-2.fa')
sarscov2_gb = os.path.join(datadir, 'Sars-Cov-2.gb')

preset_genomes = {
           'Mbovis-AF212297':{'sequence':mbovis_genome, 'gb':mbovis_gb, 'mask':mbovis_mask},
           'MTB-H37Rv':{'sequence':mtb_genome, 'gb':mtb_gb, 'mask':mtb_mask},
           'MAP-K10':{'sequence':map_genome, 'gb':map_gb,'mask':map_mask},
           'M.smegmatis-MC2155':{'sequence':msmeg_genome, 'gb':msmeg_gb},
           'Mycoplasmabovis-PG45':{'sequence':mycoplasmabovis_genome, 'gb':mycoplasmabovis_gb},
           'Sars-Cov-2':{'sequence':sarscov2_genome, 'gb':sarscov2_gb}
           }

#windows only path to binaries
bin_path = os.path.join(config_path, 'binaries')
#this is a custom filter
default_filter = 'QUAL>=40 && FORMAT/DP>=30 && DP4>=4'
annotatestr = '"AD,ADF,ADR,DP,SP,INFO/AD,INFO/ADF,INFO/ADR"'

if not os.path.exists(config_path):
    try:
        os.makedirs(config_path, exist_ok=True)
    except:
        os.makedirs(config_path)

defaults = {'threads':4, 'labelsep':'_', 'labelindex':0,
            'manifest': None,
            'unmapped':False, 'quality':25, #'trim':False
            'aligner': 'bwa', 'platform': 'illumina', 'species': None,
            'filters': default_filter, 'proximity': 10, 'mask': None,
            'uninformative_sites': False,
            'reference': None, 'gb_file': None, 'overwrite':False,
            'subset': None, 'get_stats':False,
            'old_method': False,
            'buildtree':False, 'bootstraps':100}

def check_platform():
    """See if we are running in Windows"""

    if platform.system() == 'Windows':
        print('checking binaries are present')
        fetch_binaries()
    return

def copy_ref_genomes():
    """Copy default ref genome files to config dir"""

    files = glob.glob(os.path.join(datadir, '*.fa'))
    path = sequence_path
    if not os.path.exists(path):
        os.makedirs(path,exist_ok=True)
    for src in files:
        dest = os.path.join(path, os.path.basename(src))
        shutil.copy(src, dest)
    return

copy_ref_genomes()

def fetch_binaries():
    """Get windows binaries -- windows only"""

    url = "https://github.com/dmnfarrell/snipgenie/raw/master/win_binaries/"
    os.makedirs(bin_path, exist_ok=True)
    names = ['bcftools.exe','bwa.exe','samtools.exe','tabix.exe',
             'subread-align.exe','subread-buildindex.exe','fasttree.exe',
             'makeblastdb.exe','minimap2.exe','rush.exe',
             'msys-2.0.dll','msys-bz2-1.dll','msys-lzma-5.dll','msys-ncursesw6.dll','msys-z.dll']
    for n in names:
        filename = os.path.join(bin_path,n)
        if os.path.exists(filename):
            continue
        print ('fetching %s' %n)
        link = os.path.join(url,n)
        print (filename,link)
        urllib.request.urlretrieve(link, filename)
    return

def get_files_from_paths(paths, ext='*.f*q.gz', filter_list=None):
    """Get files in multiple paths.
    Args:
        ext: wildcard for file types to parse eg. *.f*q.gz]
        filter_list: list of labels that should be present in the filenames, optional
    """

    if not type(paths) == list:
        paths = [paths]
    files=[]
    for path in paths:
        if not os.path.exists(path):
            print ('the folder %s does not exist' %path)
        s = glob.glob(os.path.join(path,'**/'+ext), recursive=True)
        files.extend(s)
    found = []
    if filter_list != None:
        for f in files:
            for n in filter_list:
                if n in f:
                    found.append(f)
        files=found
    return files

def get_samples(filenames, sep='-', index=0):
    """Get sample pairs from list of files, usually fastq. This
     returns a dataframe of unique sample labels for the input and tries
     to recognise the paired files.
     Args:
        sep: separator to split name on
        index: placement of label in split list, default 0
    Returns:
        a dataframe of samples, one row per filename
     """

    res = []
    cols = ['name','sample','filename']
    for filename in filenames:
        name = os.path.basename(filename)#.split('.')[0]
        name = name.removesuffix('.fastq.gz')
        #make sure we remove pair numbers at end before getting sample
        if name[-2:] == '_1' or name[-2:] == '_2':
            label = name[:-2]
        else:
            label = name
        #print (label)
        sample = label.split(sep)[index]
        x = [name, sample, os.path.abspath(filename)]
        res.append(x)

    df = pd.DataFrame(res, columns=cols)
    df = df.sort_values(['sample','filename'])
    df['pair'] = df.groupby('sample').cumcount()+1
    df = df.drop_duplicates('filename')
    return df

def find_best_separator(files):
    """
    Finds the best separator or substring to generate unique prefixes for a list of filenames.
    Ensures the separator produces meaningful and unique prefixes.

    Parameters:
        files (list): List of filenames (strings).

    Returns:
        str: The best separator or substring that generates unique prefixes, or None if no solution is found.
    """
    separators = ['_', '-', '.']

    # Test predefined separators
    for separator in separators:
        prefixes = [file.split(separator)[0] + separator for file in files]
        # Ensure prefixes are unique and meaningful (not the full file name)
        if len(prefixes) == len(set(prefixes)) and all(len(prefix) < len(file) for prefix, file in zip(prefixes, files)):
            return separator

    # Fallback: Test substrings for uniqueness
    for i in range(1, len(files[0])):
        # Extract substrings of length i
        prefixes = [file[:i] for file in files]
        # Check for uniqueness
        if len(prefixes) == len(set(prefixes)):
            # Return the substring as the separator
            return files[0][:i]
    return None

def get_pivoted_samples(df, allow_errors=False):
    """
    Args:
        df: dataframe of sample names from get_samples()
    Returns:
        dataframe of pivoted samples
    Get pivoted samples by pair, returns a table with one sample per row and
    filenames in separate columns. If there are duplicates they will be printed
    out and an error called. The return value is None if there is an error.
    """

    #check pairs
    c = df['sample'].value_counts()
    c = c[c>2]
    if len(c)>0:
        print ('there seem to be duplicates:')
        print (c)
        print ('here are the duplicate entries:')
        print (df[df['sample'].isin(c.index)])
        print ()
    p = pd.pivot_table(df,index='sample',columns='pair',values=['filename','name'],
                        aggfunc='first')
    if len(p.columns) > 4:
        print ('error in filename parsing, check for duplicates or labelsep/labelindex options')
        #print (p[:2])
        if allow_errors == False:
            return
    c = list(zip(p.columns.get_level_values(0),p.columns.get_level_values(1)))
    p.columns = [i[0]+str(i[1]) for i in c]
    p = p.fillna('')
    p = p.reset_index()
    return p

def check_samples_unique(samples):
    """Check that sample names are unique"""

    x = samples[samples['sample'].duplicated()]
    if len(x)>0:
        return False

def write_samples(df, path):
    """Write out sample names only using dataframe from get_samples"""

    filename = os.path.join(path, 'samples.txt')
    df.to_csv(filename, index=False, header=False)
    return filename

def get_samples_from_bams(filenames, sep='_', index=0):
    """Samples names from a list of bam files"""

    res = []
    cols = ['name','sample','bam_file']
    for filename in filenames:
        name = os.path.basename(filename)
        label = name.removesuffix('.bam').split(sep)[index]
        x = [name, label, os.path.abspath(filename)]
        res.append(x)

    df = pd.DataFrame(res, columns=cols)
    df = df.drop_duplicates('bam_file')
    return df

def check_samples_aligned(samples, outdir):
    """Check how many samples already aligned"""

    found = glob.glob(os.path.join(outdir,'*.bam'))
    print ('%s/%s samples already aligned' %(len(found),len(samples)))
    return

def mapping_stats(samples):
    """Get stats on mapping of samples"""

    def get_stats(x):
        d = tools.samtools_flagstat(x)
        s = pd.Series(d)
        return s

    for i,r in samples.iterrows():
        s = get_stats(r.bam_file)

        samples.loc[i,'mapped'] = s['primary']
        total = tools.get_fastq_read_lengths(r.filename1)
        #if 'filename2' in samples.columns:
        #    total += tools.get_fastq_read_lengths(r.filename1)
        samples.loc[i,'reads'] = total
        samples.loc[i,'perc_mapped'] = round(s['primary']/(total*2)*100,2)
    return samples

def clean_folders(samples, path):
    """Clean folders from variant_calling dir that are not in samples.
    Use with caution."""

    samplenames = list(samples['sample'])
    folders = (glob.glob(path+'/*'))
    for f in folders:
        if not os.path.isdir(f):
            continue
        name = os.path.basename(f)
        if name not in samplenames:
            print (name)
            shutil.rmtree(f)
    return

def clean_bam_files(samples, path, remove=False):
    """Check if any bams in output not in samples and remove. Not used in workflow."""

    bams = get_files_from_paths(os.path.abspath(path), '*.bam')
    print ('%s bam files and %s samples found' %(len(bams),len(samples)))
    found = set(bams)-set(samples.bam_file)
    print ('bam files no longer present in samples:')
    print (found)
    if remove == True:
        for f in found:
            print ('removed %s' %f)
            os.remove(f)
    return

def fetch_contam_file():
    """Get contam sequences"""

    url = "https://github.com/dmnfarrell/snipgenie/raw/master/extra/contam.fa.gz"
    os.makedirs(bin_path, exist_ok=True)

    tempdir = tempfile.gettempdir()
    filename = os.path.join(tempdir,'contam.fa.gz')
    destfile = os.path.join(sequence_path,'contam.fa')
    if os.path.exists(destfile):
        return
    print ('fetching contaminant sequences..')
    #link = os.path.join(url,n)
    print (filename)
    urllib.request.urlretrieve(url, filename)
    tools.gunzip(filename, destfile)
    return

def blast_contaminants(filename, limit=2000, random=False, pident=98, qcovs=90):
    """Blast reads to contaminants database
    Returns: percentages of reads assigned to each species.
    """

    fetch_contam_file()
    path = os.path.join(sequence_path,'contam.fa')
    tools.make_blast_database(path)
    if random == True:
        seqs = tools.fastq_random_seqs(filename, limit)
    else:
        seqs = tools.fastq_to_rec(filename, limit)
    bl = tools.blast_sequences(path,seqs,maxseqs=1)
    bl['stitle'] = bl.stitle.apply(lambda x: x.split('__')[0])
    bl = bl[(bl.qcovs>qcovs) & (bl.pident>pident)]
    c = bl.stitle.value_counts()
    c = pd.DataFrame(c)
    c.columns = ['hits']
    c['perc_hits'] = c.hits/limit*100
    c=c[c.hits>5]
    #print (c)
    return c

def align_reads(df, idx, outdir='mapped', callback=None, aligner='bwa', platform='illumina',
                unmapped=None, stats=True, **kwargs):
    """
    Align multiple files. Requires a dataframe with a 'sample' column to indicate
    paired files grouping. If a trimmed column is present these files will align_reads
    instead of the raw ones.
    Args:
        df: dataframe with sample names and filenames
        idx: index name
        outdir: output folder
        unmapped: folder for unmapped files if required
    """

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    if unmapped != None and not os.path.exists(unmapped):
        os.makedirs(unmapped, exist_ok=True)

    new = []
    samtoolscmd = tools.get_cmd('samtools')
    for i,r in df.iterrows():
        name = r['sample']
        if 'trimmed1' in df.columns:
            print('using trimmed files..')
            file1 = r.trimmed1
            file2 = r.trimmed2
        else:
            file1 = r.filename1
            if 'filename2' in df.columns:
                file2 = r.filename2
                #print (file2, type(file2))
            else:
                file2 = None

        out = os.path.join(outdir,str(name)+'.bam')
        #print (name)
        if aligner == 'bwa':
            aligners.bwa_align(file1, file2, idx=idx, out=out, unmapped=unmapped, **kwargs)
        elif aligner == 'bowtie':
            idx = os.path.splitext(os.path.basename(idx))[0]
            aligners.bowtie_align(file1, file2, idx=idx, out=out, **kwargs)
        elif aligner == 'subread':
            idx = os.path.splitext(os.path.basename(idx))[0]
            aligners.subread_align(file1, file2, idx=idx, out=out, **kwargs)
        elif aligner == 'minimap2':
            #idx = os.path.splitext(os.path.basename(idx))[0]
            aligners.minimap2_align(file1, file2, idx=idx, out=out, platform=platform, **kwargs)
        bamidx = out+'.bai'
        if not os.path.exists(bamidx) or kwargs['overwrite']==True:
            #print('indexing %s' %name)
            cmd = '{s} index {o}'.format(o=out,s=samtoolscmd)
            subprocess.check_output(cmd,shell=True)
            #print (cmd)
        #set bam file
        df.loc[i,'bam_file'] = os.path.abspath(out)

        #find mean depth/coverage
        if stats==True and ('meandepth' not in df.columns or pd.isnull(df.loc[i,'meandepth'])):
            cols = ['coverage','meandepth']
            c = tools.samtools_coverage(out)
            df.loc[i,cols] = c[cols]

    return df

def mpileup(bam_file, ref, out, overwrite=False):
    """Run bcftools for single file."""

    bcftoolscmd = tools.get_cmd('bcftools')
    if os.path.exists(out):
        return

    cmd = '{bc} mpileup -a {a} -O b --min-MQ 60 -o {o} -f {r} {b}'\
            .format(r=ref, b=bam_file, o=out, bc=bcftoolscmd, a=annotatestr)
    #print (cmd)
    subprocess.check_output(cmd, shell=True)
    cmd = 'bcftools index {o}'.format(o=out)
    subprocess.check_output(cmd, shell=True)
    return

def mpileup_region(region,out,bam_files,callback=None):
    """Run bcftools for single region."""

    bcftoolscmd = tools.get_cmd('bcftools')
    cmd = 'bcftools mpileup -r {reg} -O b -o {o} -f {r} {b}'.format(r=ref, reg=region, b=bam_files, o=out)
    if callback != None:
        callback(cmd)
    subprocess.check_output(cmd, shell=True)
    cmd = 'bcftools index {o}'.format(o=out)
    subprocess.check_output(cmd, shell=True)
    return

def worker(args):
    mpileup(args[0], args[1], args[2])

def mpileup_parallel_old(bam_files, ref, outpath, threads=4, callback=None, tempdir=None):
    """Run mpileup in over multiple regions with GNU parallel on linux or rush on Windows
      Separate bcf files are then joined together.
      Assumes alignment to a bacterial reference with a single chromosome.
    """

    if tempdir == None:
        tempdir = tempfile.tempdir
    bam_files = ' '.join(bam_files)
    rawbcf = os.path.join(outpath,'raw.bcf')
    chr = tools.get_chrom(ref)
    length = tools.get_fasta_length(ref)
    x = np.linspace(1,length,threads+1,dtype=int)
    print (x)

    #split genome into blocks
    blocks=[]
    for i in range(len(x)):
        if i < len(x)-1:
            blocks.append((x[i],x[i+1]-1))

    #get temp outfile names
    outfiles = []
    regions = []
    for start,end in blocks:
        region = '"{c}":{s}-{e}'.format(c=chr,s=start,e=end)
        regions.append(region)
        out = os.path.join(tempdir,'{s}-{e}.bcf'.format(s=start,e=end))
        outfiles.append(out)

    regstr = ' '.join(regions)
    #print (regstr)
    filesstr = ' '.join(outfiles)
    bcftoolscmd = tools.get_cmd('bcftools')

    if platform.system() == 'Windows':
        rushcmd = tools.get_cmd('rush')
        cmd = 'echo {reg} | {rc} -D " " "{bc} mpileup -r {{}} -f {r} -a {a} --max-depth 500 --min-MQ 60 {b} -o {p}/{{@[^:]*$}}.bcf"'\
                .format(rc=rushcmd,bc=bcftoolscmd,reg=regstr,r=ref,b=bam_files,a=annotatestr,p=tempdir)
    else:
        cmd = 'parallel bcftools mpileup -r {{1}} -a {a} -O b --max-depth 500 --min-MQ 60 -o {{2}} -f {r} {b} ::: {reg} :::+ {o}'\
                .format(r=ref, reg=regstr, b=bam_files, o=filesstr, a=annotatestr)
    print (cmd)
    #if callback != None:
    #    callback(cmd)
    subprocess.check_output(cmd, shell=True)
    #concat the separate files
    cmd = '{bc} concat {i} -O b -o {o}'.format(bc=bcftoolscmd,i=' '.join(outfiles),o=rawbcf)
    print (cmd)
    subprocess.check_output(cmd, shell=True)
    #remove temp files
    for f in outfiles:
        os.remove(f)
    return rawbcf

def mpileup_parallel(bam_file, ref, rawbcf, threads=4, callback=None, tempdir=None, show_cmd=False):
    """Run mpileup in over multiple regions with GNU parallel on linux or rush on Windows
      Separate bcf files are then joined together.
      Assumes alignment to a bacterial reference with a single chromosome.
    """

    if tempdir == None:
        tempdir = tempfile.tempdir
    #bam_files = ' '.join(bam_files)
    #rawbcf = os.path.join(outpath,'raw.bcf')
    chr = tools.get_chrom(ref)
    length = tools.get_fasta_length(ref)
    x = np.linspace(1,length,threads+1,dtype=int)
    print (x)

    #split genome into blocks
    blocks=[]
    for i in range(len(x)):
        if i < len(x)-1:
            blocks.append((x[i],x[i+1]-1))

    #get temp outfile names
    outfiles = []
    regions = []
    for start,end in blocks:
        region = '"{c}":{s}-{e}'.format(c=chr,s=start,e=end)
        regions.append(region)
        out = os.path.join(tempdir,'{s}-{e}.bcf'.format(s=start,e=end))
        outfiles.append(out)

    regstr = ' '.join(regions)
    #print (regstr)
    filesstr = ' '.join(outfiles)
    bcftoolscmd = tools.get_cmd('bcftools')
    cmd = 'parallel bcftools mpileup -r {{1}} -a {a} -O b --max-depth 500 --min-MQ 60 -o {{2}} -f {r} {b} ::: {reg} :::+ {o}'\
                .format(r=ref, reg=regstr, b=bam_file, o=filesstr, a=annotatestr)

    subprocess.check_output(cmd, shell=True, stderr=subprocess.PIPE)
    #concat the separate files
    cmd = '{bc} concat {i} --threads {t} -O b -o {o}'.format(bc=bcftoolscmd,i=' '.join(outfiles),o=rawbcf,t=threads)
    if show_cmd == True:
        print (cmd)
    subprocess.check_output(cmd, shell=True, stderr=subprocess.PIPE)
    #remove temp files
    for f in outfiles:
        os.remove(f)
    return rawbcf

def variant_calling_old(bam_files, ref, outpath, samples=None, threads=4,
                    callback=None, overwrite=False, filters=None, gff_file=None,
                    mask=None, tempdir=None, sep='_',
                    proximity=10, **kwargs):
    """Call variants with bcftools"""

    st = time.time()
    if samples is not None:
        #write out new samples.txt file for reheader step
        write_samples(samples[['sample']], outpath)
        sample_file = os.path.join(outpath,'samples.txt')
    if filters == None:
        filters = default_filter
    rawbcf = os.path.join(outpath,'raw.bcf')
    bcftoolscmd = tools.get_cmd('bcftools')
    if not os.path.exists(rawbcf) or overwrite == True:
        print ('running mpileup for %s files..' %len(bam_files))
        if threads == 1:
            bam_files = ' '.join(bam_files)
            cmd = '{bc} mpileup -a {a} --max-depth 500 -O b --min-MQ 60 -o {o} -f {r} {b}'\
                .format(bc=bcftoolscmd,r=ref, b=bam_files, o=rawbcf, a=annotatestr)
            #print (cmd)
            tmp = subprocess.check_output(cmd, shell=True)
        #or use mpileup in parallel to speed up
        else:
            rawbcf = mpileup_parallel_old(bam_files, ref, outpath, threads=threads,
                                        tempdir=tempdir, callback=callback)
        #new method
        #rawbcf = mpileup_multiprocess(bam_files, ref, outpath, threads=threads,
        #                                 callback=callback)

    else:
        print ('%s already exists' %rawbcf)
        #check existing file samples here
        rawsamples = tools.get_vcf_samples(rawbcf)
        samples = pd.read_csv(sample_file,names=['name'])
        if len(samples) != len(rawsamples):
            print ('WARNING: samples in raw.bcf appear to be different to current samples.')
            print ('You may have added files since the previous run and will need to overwrite raw.bcf')

    #find snps only
    print ('calling variants..')
    vcfout = os.path.join(outpath,'calls.vcf')
    cmd = '{bc} call --ploidy 1 -m -v -o {o} {raw}'.format(bc=bcftoolscmd,o=vcfout,raw=rawbcf)
    print (cmd)
    subprocess.check_output(cmd,shell=True)
    sites = tools.bcftools_count_sites(vcfout)
    print ('%s sites called as variants' %sites)

    #relabel samples in vcf header
    if samples is not None:
        relabel_vcfheader(vcfout, sample_file)

    #filters
    filtered = os.path.join(outpath,'filtered.vcf.gz')
    cmd = '{bc} filter -i "{f}" -o {o} -O z {i}'.format(bc=bcftoolscmd,i=vcfout,o=filtered,f=filters)
    print (cmd)
    tmp = subprocess.check_output(cmd,shell=True)

    #get only snps
    print ('splitting snps and indels..')
    snpsout = os.path.join(outpath,'snps.vcf.gz')
    cmd = '{bc} view -v snps -o {o} -O z {i}'.format(bc=bcftoolscmd,o=snpsout,i=filtered)
    print (cmd)
    subprocess.check_output(cmd,shell=True)

    #also get indels only to separate file
    indelsout = os.path.join(outpath,'indels.vcf.gz')
    #cmd = '{bc} call -V snps --ploidy 1 -m -v -o {o} {raw}'.format(bc=bcftoolscmd,o=indelsout,raw=rawbcf)
    cmd = '{bc} view -v indels -o {o} -O z {i}'.format(bc=bcftoolscmd,o=indelsout,i=filtered)
    print (cmd)
    subprocess.check_output(cmd,shell=True)

    #apply mask if provided
    if mask != None:
        mask_filter(snpsout, mask, outdir=outpath, overwrite=True)
        #mask_filter(indelsout, mask, outdir=outpath, overwrite=True)

    #prox filter
    site_proximity_filter(snpsout, dist=int(proximity), outdir=outpath, overwrite=True)
    #consequence_calling
    consequence_calling(snpsout, indelsout, outpath, gff_file, ref)
    print ('took %s seconds' %str(round(time.time()-st,0)))
    return snpsout

def consequence_calling(snpsout, indelsout, outpath, gff_file, ref):
    """Consequence calling using gff file"""

    if gff_file != None:
        print ('consequence calling..')
        try:
            csqout = os.path.join(outpath, 'csq.tsv')
            m = csq_call(ref, gff_file, snpsout, csqout)
            m.to_csv(os.path.join(outpath,'csq.matrix'))
            #indels as well
            if indelsout is not None:
                csqout = os.path.join(outpath, 'csq_indels.tsv')
                m = csq_call(ref, gff_file, indelsout, csqout)
                m.to_csv(os.path.join(outpath,'csq_indels.matrix'))
        except Exception as e:
            print (e)
    return

#new code

def run_file(bam_file, ref, outpath, overwrite=False, threads=4):
    """
    Run mpileup/snp calling for one sample.
    """

    if not os.path.exists(outpath):
        os.makedirs(outpath, exist_ok=True)

    bcftoolscmd = tools.get_cmd('bcftools')
    rawbcf = os.path.join(outpath,'raw.bcf')

    if not os.path.exists(rawbcf) or overwrite == True:
        mpileup_parallel(bam_file, ref, rawbcf, threads=threads)

    #call
    vcfout = os.path.join(outpath,'calls.vcf')
    snpsout = os.path.join(outpath,'snps.bcf')
    indelsout = os.path.join(outpath,'indels.bcf')

    #assumes other files are there and returns - improve this
    if os.path.exists(vcfout) and overwrite == False:
        return snpsout, indelsout
    tools.bcftools_call(rawbcf, vcfout, show_cmd=False)

    #print ('splitting snps and indels..')
    cmd = '{bc} view -v snps -o {o} -O b {i}'.format(bc=bcftoolscmd,o=snpsout,i=vcfout)
    #print (cmd)
    subprocess.check_output(cmd,shell=True)
    cmd = 'bcftools index {o}'.format(o=snpsout)
    subprocess.check_output(cmd, shell=True)

    #also get indels only to separate file
    cmd = '{bc} view -v indels -o {o} -O b {i}'.format(bc=bcftoolscmd,o=indelsout,i=vcfout)
    #print (cmd)
    subprocess.check_output(cmd,shell=True)
    cmd = 'bcftools index {o}'.format(o=indelsout)
    subprocess.check_output(cmd, shell=True)
    return snpsout, indelsout

def write_filelist(files, filename):
    """Write out sample names only using dataframe from get_samples"""

    df = pd.DataFrame(files)
    df.to_csv(filename, index=False, header=False)
    return

def smart_tqdm(iterable, **kwargs):
    """check if we are in notebook"""
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            from tqdm.notebook import tqdm
        else:
            from tqdm import tqdm
    except NameError:
        from tqdm import tqdm
    return tqdm(iterable, **kwargs)

def variant_calling(samples, ref, outpath, filters=None, proximity=10, mask=None,
                    gff_file=None, overwrite=False, threads=4, **kwargs):
    """
    Run variant calling for multiple samples.
    """

    bcftoolscmd = tools.get_cmd('bcftools')
    vcfoutpath = os.path.join(outpath, 'variant_calling')
    os.makedirs(vcfoutpath, exist_ok=True)
    if samples is not None:
        #write out new samples.txt file for reheader step
        write_samples(samples[['sample']], outpath)
        sample_file = os.path.join(outpath,'samples.txt')

    outfiles = []
    indelfiles = []
    #from tqdm import tqdm
    for i, r in smart_tqdm(samples.iterrows(), total=len(samples)):
    #for i,r in samples.iterrows():
        name = r['sample']
        #print (name)
        path = os.path.join(vcfoutpath, name)
        #print (path, os.path.exists(path))
        snpfile,indfile = run_file(r.bam_file, ref, path, threads=threads)
        outfiles.append(snpfile)
        indelfiles.append(indfile)

    #write file list
    bcffilelist = os.path.join(outpath, 'bcffiles.txt')
    write_filelist(outfiles, bcffilelist)
    merged = os.path.join(outpath, 'merged.vcf.gz')
    #bcf_files = ' '.join(outfiles)
    #print (outfiles)
    cmd = '{bc} merge --threads {t} -0 -O z -o {o} --file-list {f}'.format(f=bcffilelist,
                                                                           o=merged, bc=bcftoolscmd,t=threads)
    print (cmd)
    subprocess.check_output(cmd, shell=True)
    relabel_vcfheader(merged, sample_file)

    #filter snp calls
    if filters == None:
        filters = default_filter
    if filters != '':
        filtered = os.path.join(outpath,'snps.vcf.gz')
        cmd = '{bc} filter -i "{f}" -o {o} -O z {i}'.format(bc=bcftoolscmd,i=merged,
                                                            o=filtered,f=filters)
        print (cmd)
        tmp = subprocess.check_output(cmd,shell=True)
    else:
        filtered = merged

    #apply mask if required
    if mask != None:
        mask_filter(filtered, mask, outdir=outpath, overwrite=True)
        #mask_filter(indelsout, mask, outdir=outpath, overwrite=True)
    #proximity filter
    site_proximity_filter(filtered, dist=int(proximity), outdir=outpath, overwrite=True)

    #handle indels
    indelfilelist = os.path.join(outpath, 'indelfiles.txt')
    write_filelist(indelfiles, indelfilelist)
    indelsout = os.path.join(outpath, 'indels.vcf.gz')
    cmd = '{bc} merge --threads {t} -0 -O z -o {o} --file-list {f}'.format(f=indelfilelist,
                                                                           o=indelsout, bc=bcftoolscmd,t=threads)
    print (cmd)
    tmp = subprocess.check_output(cmd,shell=True)
    relabel_vcfheader(indelsout, sample_file)
    #consequence_calling
    consequence_calling(filtered, indelsout, outpath, gff_file, ref)
    #print ('took %s seconds' %str(round(time.time()-st,0)))
    return filtered

def csq_call(ref, gff_file, vcf_file, csqout):
    """Consequence calling"""

    bcftoolscmd = tools.get_cmd('bcftools')
    cmd = '{bc} csq -f {r} -g {g} {f} -v 0 -Ot -o {o}'.format(bc=bcftoolscmd,r=ref,g=gff_file,
                f=vcf_file,o=csqout)
    print (cmd)
    #if callback != None:
    #    callback(cmd)
    tmp = subprocess.check_output(cmd,shell=True)
    csqdf = read_csq_file(csqout)
    #get presence/absence matrix of csq mutations
    m = get_aa_snp_matrix(csqdf)
    return m

def relabel_vcfheader(vcf_file, sample_file):
    """Re-label samples in vcf header"""

    bcftoolscmd = tools.get_cmd('bcftools')
    rlout = os.path.join(tempdir,'calls.vcf')
    cmd = '{bc} reheader --samples {s} -o {o} {v}'.format(bc=bcftoolscmd,o=rlout,
                                                v=vcf_file,s=sample_file)
    print(cmd)
    tmp = subprocess.check_output(cmd,shell=True)
    #rewrite file
    shutil.copy(rlout, vcf_file)
    #remove temp file
    os.remove(rlout)
    return

def mask_filter(vcf_file, mask_file, overwrite=False, outdir=None):
    """Remove any masked sites using a bed file, overwrites input"""

    print('using mask bed file', mask_file)
    mask = pd.read_csv(mask_file,sep='\t',names=['chrom','start','end'])
    #print (mask)
    def do_mask(x,i):
        #print (x.start, x.end, i)
        if (x.start<=i) & (x.end>=i):
            return 1
    import vcf
    print (vcf_file)
    vcf_reader = vcf.Reader(open(vcf_file, 'rb'))
    sites = [record.POS for record in vcf_reader]
    found = []
    for i in sites:
        m = mask.apply( lambda x: do_mask(x,i),1)
        m = m[m==1]
        if len(m)>0:
            #print (i)
            found.append(i)
    print('found %s/%s sites in masked regions' %(len(found),len(sites)))
    new = sorted(list(set(sites) - set(found)))
    if outdir == None:
        outdir = tempfile.gettempdir()
    if overwrite == True:
        overwrite_vcf(vcf_file, new, outdir)
    return

def site_proximity_filter(vcf_file, dist=10, overwrite=False, outdir=None):
    """Remove any pairs of sites within dist of each other.
    Args:
        vcf_file: input vcf file with positions to filter
        dist: distance threshold
        overwrite: whether to overwrite the vcf
    """

    if dist == 0:
        return
    print('applying proximity filter..')
    #get vcf positions into dataframe
    df = tools.get_vcf_positions(vcf_file)
    #print (df)
    sites = list(df.pos.unique())
    found = []
    #check distances in sites per sample
    for s, g in df.groupby(['sample']):
        pos = list(g.pos)
        for i in range(len(pos)-1):
            if pos[i+1] - pos[i] <= dist:
                found.extend([pos[i], pos[i+1]])
    #all unique positions
    found = list(set(found))
    new = sorted(list(set(sites) - set(found)))
    print ('proximity filter removed %s/%s sites' %(len(found),len(sites)))
    if overwrite == True:
        overwrite_vcf(vcf_file, new, outdir)
    return

def overwrite_vcf(vcf_file, sites, outdir=None):
    """Make a new vcf with subset of sites"""

    if outdir == None:
        outdir = tempfile.gettempdir()

    import vcf
    out = os.path.join(outdir,'temp.vcf')
    vcf_reader = vcf.Reader(open(vcf_file, 'rb'))
    vcf_writer = vcf.Writer(open(out, 'w'), vcf_reader)
    for record in vcf_reader:
        if record.POS in sites:
            #print (record)
            vcf_writer.write_record(record)
    vcf_writer.close()
    #copy or overwrite input vcf
    bcftoolscmd = tools.get_cmd('bcftools')
    cmd = '{b} view {o} -O z -o {gz}'.format(b=bcftoolscmd,o=out,gz=vcf_file)
    print (cmd)
    tmp = subprocess.check_output(cmd,shell=True)
    return

def trim_files(df, outpath, overwrite=False, threads=4, quality=30):
    """Batch trim fastq files"""

    method = 'cutadapt'
    if platform.system() == 'Windows':
        method = 'default'
    if not os.path.exists(outpath):
        os.makedirs(outpath, exist_ok=True)
    for i,row in df.iterrows():
        out1 = os.path.join(outpath, os.path.basename(row.filename1))
        out2 = os.path.join(outpath, os.path.basename(row.filename2))
        if not os.path.exists(out1) or overwrite == True:
            out1,out2 = tools.trim_reads(row.filename1, row.filename2,
                    outpath, threads=threads, quality=quality, method=method)

        df.loc[i,'trimmed1'] = out1
        df.loc[i,'trimmed2'] = out2
    return df

def read_csq_file(filename):
    """Read csq tsv outpt file into dataframe"""

    cols = ['1','sample','2','chrom','start','snp_type','gene','locus_tag','strand','feature_type','aa','nuc',]
    csqdf = pd.read_csv(filename,sep='[|\t]',comment='#',names=cols, engine='python')
    csqdf = csqdf[csqdf['1']!='LOG']
    csqdf['aa'] = csqdf.aa.fillna(csqdf.snp_type)
    csqdf['nuc'] = csqdf.nuc.fillna(csqdf.snp_type)
    return csqdf

def get_aa_snp_matrix(df):
    """Get presence/absence matrix from csq calls table"""

    df = df.drop_duplicates(['gene','aa','sample'])
    x = df.set_index(['start','gene','aa','snp_type','sample'])['nuc'].unstack('sample')
    x[x.notna()] = 1
    x = x.fillna(0)
    return x

def run_bamfiles(bam_files, ref, gff_file=None, mask=None, outdir='.', threads=4,
                    sep='_', labelindex=0, old_method=False, **kwargs):
    """
    Run workflow with bam files from a previous sets of alignments.
    We can arbitrarily combine results from multiple other runs this way.
    kwargs are passed to variant_calling method.
    Should write a samples.txt file in the outdir if vcf header is to be
    relabelled.
    Args:
        bam_files: list of bam files
        ref: reference genome
        samples: dataframe of sample names, if not provided try to get from bam files
        sep: separator for getting sample names, will write out a sample.txt file
         in the output folder for bcftools to usef for relabelling the vcf
    """

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    print ('%s samples were loaded:' %len(bam_files))

    samples = get_samples_from_bams(bam_files, sep=sep)
    if old_method == True:
        vcf_file = variant_calling_old(bam_files, ref, outdir,
                                    samples=samples, threads=threads,
                                    gff_file=gff_file, mask=mask, sep=sep,
                                    **kwargs)
    else:
        vcf_file = variant_calling(samples, ref, outdir,
                                    threads=threads,
                                    gff_file=gff_file, mask=mask, sep=sep,
                                    **kwargs)
    run_vcf(vcf_file, outdir, threads=threads)
    return

def run_vcf(vcf_file, outdir, threads=8):
    """Run steps on a final vcf file to get core alignment, snp matrix and tree"""

    print ('getting core alignment..')
    snprecs, smat = tools.core_alignment_from_vcf(vcf_file)
    outfasta = os.path.join(outdir, 'core.fa')
    SeqIO.write(snprecs, outfasta, 'fasta')
    smat.to_csv(os.path.join(outdir,'core.txt'), sep=' ')
    aln = AlignIO.read(outfasta, 'fasta')
    #remove ref
    aln = aln[1:]
    print ('computing snp distance matrix..')
    snp_dist = tools.snp_dist_matrix(aln)
    snp_dist.to_csv(os.path.join(outdir,'snpdist.csv'), sep=',')
    #treefile = trees.run_RAXML(outfasta, outpath=outdir, threads=threads)
    #ls = len(smat)
    #trees.convert_branch_lengths(treefile,os.path.join(outdir,'tree.newick'), ls)
    return

class Logger(object):
    """
    This class duplicates sys.stdout to a log file
    source: https://stackoverflow.com/q/616645
    """
    def __init__(self, filename="run.log", mode="a"):
        self.stdout = sys.stdout
        self.file = open(filename, mode)
        sys.stdout = self

    def __del__(self):
        self.close()

    def __enter__(self):
        pass

    def __exit__(self, *args):
        self.close()

    def write(self, message):
        self.stdout.write(message)
        self.file.write(message)

    def flush(self):
        self.stdout.flush()
        self.file.flush()
        os.fsync(self.file.fileno())

    def close(self):
        if self.stdout != None:
            sys.stdout = self.stdout
            self.stdout = None

        if self.file != None:
            self.file.close()
            self.file = None

class WorkFlow(object):
    """Class for implementing a prediction workflow from a set of options"""
    def __init__(self, **kwargs):
        for i in kwargs:
            self.__dict__[i] = kwargs[i]
        for i in defaults:
            if i not in self.__dict__:
                self.__dict__[i] = defaults[i]
        #make output folder
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir, exist_ok=True)
        #start logger
        self.logfile = os.path.join(self.outdir, 'run.log')
        #sys.stdout = Logger(self.logfile)
        print ('The following options were supplied')
        dt_string = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        print("time: ", dt_string)
        print ('-------')
        for i in self.__dict__:
            print (i, ':', self.__dict__[i])
        print ()
        return

    def setup(self):
        """Setup main parameters"""

        if self.species != None:
            s = self.species
            if s not in preset_genomes:
                valid = '; '.join(list(preset_genomes.keys()))
                print ('Invalid preset species value! Use one of: %s' %valid)
                return
            self.reference = preset_genomes[s]['sequence']
            self.gb_file = preset_genomes[s]['gb']
            if s == 'Mbovis-AF212297':
                self.mask = mbovis_mask
        elif self.reference == None:
            self.reference = mbovis_genome
            self.gb_file = mbovis_gb

        if self.threads == None:
            import multiprocessing
            self.threads = multiprocessing.cpu_count()
        else:
            self.threads = int(self.threads)

        if self.manifest != None:
            print ('using manifest file for samples')
            df = pd.read_csv(self.manifest)
            #ensure sample names are strings
            df['sample'] = df['sample'].astype(str)
        else:
            self.filenames = get_files_from_paths(self.input)
            #get files from input folder(s)
            if self.labelsep == 'guess':
                #guess separator
                sep = find_best_separator(self.filenames)
                if sep is not None:
                    self.labelsep = sep
                    print (f'using {sep} as separator')
                else:
                    print ('cannot guess separator, you need to supply one')
                    return

            df = get_samples(self.filenames, sep=self.labelsep, index=self.labelindex)
            df = get_pivoted_samples(df)

        if df is None:
            return
        if len(df) == 0:
            print ('no samples provided. files should be fastq.gz type')
            return False

        sample_size = len(df['sample'].unique())
        print ('%s samples were loaded:' %sample_size)
        print ('----------------------')
        print (df)
        print ()
        s = check_samples_unique(df)
        if self.get_stats == True:
            print ('getting read lengths..')
            df['read_length'] = df.filename1.apply(tools.get_fastq_info)
        self.fastq_table = df
        if s == False:
            print ('samples names are not unique! try a different labelsep value.')
            return False
        #save sample table
        df.to_csv(os.path.join(self.outdir,'samples.csv'),index=False)

        print ('building index')
        if self.aligner == 'bwa':
            aligners.build_bwa_index(self.reference)
        elif self.aligner == 'bowtie':
            aligners.build_bowtie_index(self.reference)
        elif self.aligner == 'subread':
            aligners.build_subread_index(self.reference)
        if self.gb_file != None:
            #convert annotation to gff for consequence calling
            self.gff_file = os.path.join(self.outdir, os.path.basename(self.gb_file)+'.gff')
            tools.gff_bcftools_format(self.gb_file, self.gff_file)
        else:
            self.gff_file = None
        #set temp dir
        self.tempdir = os.path.join(self.outdir, 'tmp')
        if not os.path.exists(self.tempdir):
            os.makedirs(self.tempdir)
        time.sleep(1)
        return True

    def run(self):
        """Run workflow"""

        #this master table tracks our outputs
        samples = self.fastq_table
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir, exist_ok=True)
        #write_samples(samples[['sample']], self.outdir)
        if len(samples)==0:
            print ('no samples found')
            return

        '''if self.trim == True:
            print ('trimming fastq files')
            print ('--------------------')
            trimmed_path = os.path.join(self.outdir, 'trimmed')
            samples = trim_files(samples, trimmed_path, self.overwrite,
                                  quality=self.quality, threads=self.threads)
            print ()'''
        print ('aligning files')
        print ('--------------')
        print ('Using reference genome: %s' %self.reference)
        path = os.path.join(self.outdir, 'mapped')
        if self.unmapped == True:
            unmapped = os.path.join(self.outdir, 'unmapped')
        else:
            unmapped = None
        check_samples_aligned(samples, path)
        samples = align_reads(samples, idx=self.reference, outdir=path,
                        aligner=self.aligner, platform=self.platform,
                        unmapped=unmapped, stats=self.get_stats,
                        threads=self.threads, overwrite=self.overwrite)

        #lowdepth = samples[samples.meandepth<15]
        #if len(lowdepth)>0:
        #    print ('%s samples have mean depth <15' %len(lowdepth))

        #mapping stats
        if self.get_stats == True:
            print ('getting mapping stats..')
            samples = mapping_stats(samples)
        #save sample table agaim
        samples.to_csv(os.path.join(self.outdir,'samples.csv'),index=False)

        print ()
        print ('calling variants')
        print ('----------------')
        if self.old_method == True:
            #old variant calling method
            print ('using old calling method')
            bam_files = list(samples.bam_file)
            self.vcf_file = variant_calling_old(bam_files, self.reference, self.outdir,
                                            samples=samples,
                                            threads=self.threads,
                                            gff_file=self.gff_file,
                                            filters=self.filters,
                                            mask=self.mask,
                                            proximity=self.proximity,
                                            overwrite=self.overwrite,
                                            tempdir=self.tempdir)
        else:
            self.vcf_file = variant_calling(samples, self.reference, self.outdir,
                                            threads=self.threads,
                                            gff_file=self.gff_file,
                                            filters=self.filters,
                                            mask=self.mask,
                                            proximity=self.proximity,
                                            overwrite=self.overwrite,
                                            tempdir=self.tempdir)

        print (self.vcf_file)
        print ()
        print ('making SNP matrix')
        print ('-----------------')
        snprecs, smat = tools.core_alignment_from_vcf(self.vcf_file, self.uninformative_sites)
        outfasta = os.path.join(self.outdir, 'core.fa')
        SeqIO.write(snprecs, outfasta, 'fasta')
        #write out sites matrix as txt file
        smat.to_csv(os.path.join(self.outdir,'core.txt'), sep=' ')
        print ()
        #write out pairwise snp distances
        aln = AlignIO.read(outfasta, 'fasta')
        aln = aln[1:]
        snp_dist = tools.snp_dist_matrix(aln)
        snp_dist.to_csv(os.path.join(self.outdir,'snpdist.csv'), sep=',')

        print ('Done. Sample summary:')
        print ('---------------------')
        #pd.set_option('display.max_rows', 150)
        #print (samples.drop(columns=list(samples.filter(regex='filename'))))
        print ('%s samples processed' %len(samples))
        print ()

        if self.buildtree == True:
            print ('building tree')
            print ('-------------')
            if len(samples) <= 2:
                print ('Cannot build tree, too few samples.')
                return
            if platform.system() == 'Windows':
                outfile = os.path.join(self.outdir, 'tree.newick')
                treefile = trees.run_fasttree(outfasta, outfile)
            else:
                treefile = trees.run_RAXML(outfasta, threads=self.threads,
                            bootstraps=self.bootstraps, outpath=self.outdir)
            if treefile == None:
                return
            ls = len(smat)
            trees.convert_branch_lengths(treefile,os.path.join(self.outdir,'tree.newick'), ls)
        print ()
        #check unmapped reads

        return

def test_run(path='test_run', sample_size='small', simulate=True,
             overwrite=True, aligner='bwa', out='test_run_results', threads=4):
    """Test run.
    Args:
        path: destination path for simulated data
        simulate: simulate data even if already present
        overwrite: overwrite previous outputs - re-does alignments
    """

    from . import simulate
    ref = sarscov2_genome
    newick = os.path.join(module_path, 'testing', 'sim.newick')
    if sample_size == 'small':
        fasta_file = os.path.join(module_path,'testing','phastsim_sarscov2_small.fasta')
    elif sample_size == 'medium':
        fasta_file = os.path.join(module_path,'testing','phastsim_sarscov2_medium.fasta')
    elif sample_size == 'large':
        fasta_file = os.path.join(module_path,'testing','phastsim_sarscov2_large.fasta')
    if not os.path.exists(fasta_file):
        return
    print (f'found {fasta_file}')
    #clear output folder
    if os.path.exists(out):
        shutil.rmtree(out)
    if simulate == True or not os.path.exists(path):
        os.makedirs(path,exist_ok=True)
        print ('simulating fastq files..')
        simulate.simulate_paired_end_reads(fasta_file, 150, 1e5, path)
    start = time.time()

    args = {'threads':threads, 'outdir': out, 'input': path,
            'species':'Sars-Cov-2',
            'aligner':aligner, 'filters':'QUAL>=40 && DP4>=4',
            'reference': None, 'overwrite': overwrite}
    W = WorkFlow(**args)
    st = W.setup()
    if st == True:
        W.run()
    samples = W.fastq_table
    meta = simulate.create_meta_data(samples['sample'])
    meta.to_csv(os.path.join(out,'metadata.csv'))
    print (f'test run results in {out}. you can also open the results folder in the GUI.')
    tt = time.time()-start
    print (f'took {tt} seconds')
    return tt

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='snipgenie CLI tool. https://github.com/dmnfarrell/snipgenie')
    parser.add_argument("-i", "--input", action='append', dest="input", default=[],
                        help="input folder(s)", metavar="FILE")
    parser.add_argument("-M", "--manifest", dest="manifest", default=None,
                        help="manifest file with samples, optional - overrides input", metavar="FILE")
    parser.add_argument("-r", "--reference", dest="reference", default=None,
                        help="reference genome filename", metavar="FILE")
    parser.add_argument("-S", "--species", dest="species", default=None,
                        help="set the species reference genome, overrides -r.\
                         possible values are %s" %', '.join(preset_genomes.keys()))
    parser.add_argument("-g", "--genbank_file", dest="gb_file", default=None,
                        help="annotation file, optional", metavar="FILE")
    parser.add_argument("-t", "--threads", dest="threads", default=4,
                        help="cpu threads to use")
    parser.add_argument("-sep", "--labelsep", dest="labelsep", default='_',
                        help="symbol to split the sample labels on if parsing filenames")
    parser.add_argument("-x", "--labelindex", dest="labelindex", default=0,
                        help="position to extract label in split filenames")
    parser.add_argument("-w", "--overwrite", dest="overwrite", action="store_true", default=False,
                        help="overwrite intermediate files")
    #parser.add_argument("-T", "--trim", dest="trim", action="store_true", default=False,
    #                    help="whether to trim fastq files" )
    parser.add_argument("-U", "--unmapped", dest="unmapped", action="store_true", default=False,
                        help="whether to save unmapped reads" )
    parser.add_argument("-Q", "--quality", dest="quality", default=25,
                        help="right trim quality, default 25")
    parser.add_argument("-f", "--filters", dest="filters", default=default_filter,
                        help="variant calling post-filters" )
    parser.add_argument("-m", "--mask", dest="mask", default=None,
                        help="supply mask regions from a bed file" )
    parser.add_argument("-pf", "--proximity", dest="proximity", default=10,
                        help="proximity filter value, set 0 to not apply filter")
    parser.add_argument("-u", "--uninformative", dest="uninformative_sites", action="store_true", default=False,
                        help="keep uninformative sites when calling variants" )
    parser.add_argument("-p", "--platform", dest="platform", default='illumina',
                        help="sequencing platform, change to ont if using oxford nanopore")
    parser.add_argument("-a", "--aligner", dest="aligner", default='bwa',
                        help="aligner to use, bwa, subread, bowtie or minimap2")
    parser.add_argument("-b", "--buildtree", dest="buildtree", action="store_true", default=False,
                        help="whether to build a phylogenetic tree, requires RaXML" )
    parser.add_argument("-N", "--bootstraps", dest="bootstraps", default=100,
                        help="number of bootstraps to build tree")
    parser.add_argument("-o", "--outdir", dest="outdir",
                        help="Results folder", metavar="FILE")
    parser.add_argument("-old", "--old_method", dest="old_method", action="store_true",
                        help="use old calling method")
    parser.add_argument("-q", "--qc", dest="qc", action="store_true",
                        help="QC report")
    parser.add_argument("-s", "--stats", dest="get_stats", action="store_true", default=False,
                        help="Calculate read length and mapping stats")
    parser.add_argument("-dummy", "--dummy", dest="dummy",  action="store_true",
                        default=False, help="Check samples but don't execute run")
    parser.add_argument("-T", "--test", dest="test",  action="store_true",
                        default=False, help="Test run")
    parser.add_argument("-v", "--version", dest="version", action="store_true",
                        help="Get version")

    args = vars(parser.parse_args())
    check_platform()

    if args['test'] == True:
        test_run()
    elif args['version'] == True:
        from . import __version__
        print ('snipgenie version %s' %__version__)
        print ('https://github.com/dmnfarrell/snipgenie')
    elif args['qc'] == True:
        print ('Running qc report')
        qcfile = 'qc_report.pdf'
        filenames = get_files_from_paths(args['input'])
        tools.pdf_qc_reports(filenames, qcfile)
    elif args['outdir'] == None:
        print ('No input or output folders provided. These are required.')
        print ('Example:')
        print ('snipgenie -r <reference> -i <input folder with fastq.gz files> -o <output folder>')
        print ('Use -h for more help on options.')
    else:
        W = WorkFlow(**args)
        st = W.setup()
        if st == True and args['dummy'] == False:
            W.run()

if __name__ == '__main__':
    main()
