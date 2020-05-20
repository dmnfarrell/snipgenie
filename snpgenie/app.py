#!/usr/bin/env python

"""
    snpgenie methods for cmd line tool.
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
import time, datetime
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
from . import tools, aligners, trees
import multiprocessing as mp

tempdir = tempfile.gettempdir()
home = os.path.expanduser("~")
config_path = os.path.join(home,'.config','snpgenie')
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'data')
sequence_path = os.path.join(config_path, 'genome')
annotation_path = os.path.join(config_path, 'annotation')
mbovis_genome = os.path.join(sequence_path, 'Mbovis_AF212297.fa')
mtb_genome = os.path.join(sequence_path, 'MTB-H37Rv.fa')
mbovis_gff = os.path.join(datadir, 'Mbovis_csq_format.gff')
mtb_gff = None
#windows only path to binaries
bin_path = os.path.join(config_path, 'binaries')
default_filter = 'QUAL>=40 && INFO/DP>=30 && DP4>=4 && MQ>35'

if not os.path.exists(config_path):
    try:
        os.makedirs(config_path, exist_ok=True)
    except:
        os.makedirs(config_path)

defaults = {'threads':None, 'labelsep':'_','trim':False, 'quality':25, 'filters': default_filter,
            'reference': None, 'gff_file': None, 'overwrite':False, 'buildtree':False}

def check_platform():
    """See if we are running in Windows"""

    if platform.system() == 'Windows':
        print('checking binaries are present')
        fetch_binaries()
    return

def copy_ref_genomes():
    """Copy default ref genome files to config dir"""

    files =  glob.glob(os.path.join(datadir, '*.fa'))
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

    url = "https://github.com/dmnfarrell/snpgenie/raw/master/win_binaries/"
    os.makedirs(bin_path, exist_ok=True)
    names = ['bcftools.exe','bwa.exe','samtools.exe','tabix.exe',
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

def get_files_from_paths(paths):
    """Get files in multiple paths"""

    if not type(paths) == list:
        paths = [paths]
    files=[]
    for path in paths:
        s = glob.glob(os.path.join(path,'**/*.f*q.gz'), recursive=True)
        files.extend(s)
    return files

def get_samples(filenames, sep='-'):
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

def check_samples_unique(samples):
    """Check that sample names are unique"""

    x = samples['sample'].value_counts()
    if len(x[x>2]) > 0:
        return False

def results_summary(df):
    return df.groupby('sample').first()[['name','bam_file','read_length','perc_mapped']].reset_index()

def write_samples(df, path):
    filename = os.path.join(path, 'samples.txt')
    df.drop_duplicates('sample')['sample'].to_csv(filename,index=False,header=False)

def align_reads(samples, idx, outdir='mapped', callback=None, **kwargs):
    """
    Align multiple files. Requires a dataframe with a 'sample' column to indicate
    paired files grouping. If a trimmed column is present these files will align_reads
    instead of the raw ones.
    """

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    new = []
    samtoolscmd = tools.get_cmd('samtools')
    for name,df in samples.groupby('sample'):
        #print (name)
        if callback != None:
            callback('aligning %s' %name)
        if 'trimmed' in df.columns:
            files = list(df.trimmed)
            if callback != None:
                callback('using trimmed')
        else:
            files = list(df.filename)
        if len(files) == 1:
            #unpaired reads
            files.append(None)
        out = os.path.join(outdir,name+'.bam')
        aligners.bwa_align(files[0],files[1], idx=idx, out=out, **kwargs)
        bamidx = out+'.bai'
        if not os.path.exists(bamidx) or kwargs['overwrite']==True:
            cmd = '{s} index {o}'.format(o=out,s=samtoolscmd)
            subprocess.check_output(cmd,shell=True)
            print (cmd)
        index = df.index
        samples.loc[index,'bam_file'] = out
        #find mean depth
        #depth = tools.get_bam_depth(out)
        #samples.loc[index,'depth'] = depth
        #get mapping info and add to samples samples table
        stat = tools.samtools_flagstats(out)
        perc = stat['mapped']/stat['total']*100
        samples.loc[index,'perc_mapped'] = perc
        if callback != None:
            callback(out)
    return samples

def mpileup_region(region,out,bam_files,callback=None):
    """Run bcftools for single region."""

    bcftoolscmd = 'bcftools'
    if getattr(sys, 'frozen', False):
        bcftoolscmd = tools.resource_path('bin/bcftools.exe')
    elif platform.system() == 'Windows':
        fetch_binaries()
        cmd = os.path.join('bin_path','bcftools.exe')

    cmd = 'bcftools mpileup -r {reg} -O b -o {o} -f {r} {b}'.format(r=ref, reg=region, b=bam_files, o=out)
    if callback != None:
        callback(cmd)
    subprocess.check_output(cmd, shell=True)
    cmd = 'bcftools index {o}'.format(o=out)
    subprocess.check_output(cmd, shell=True)
    return

def mpileup_multiprocess(bam_files, ref, outpath, threads=4, callback=None):
    """Run mpileup in parallel over multiple regions, then concat vcf files.
    Assumes alignment to a bacterial reference with a single chromosome."""

    bam_files = ' '.join(bam_files)
    rawbcf = os.path.join(outpath,'raw.bcf')
    tmpdir = '/tmp'
    chr = tools.get_chrom(ref)
    length = tools.get_fasta_length(ref)

    #find regions
    bsize = int(length/(threads-1))
    x = np.linspace(1,length,threads,dtype=int)
    blocks=[]
    for i in range(len(x)):
        if i < len(x)-1:
            blocks.append((x[i],x[i+1]-1))
    #print (blocks, bsize)

    pool = mp.Pool(threads)
    outfiles = []
    st = time.time()

    for start,end in blocks:
        print (start, end)
        region = '{c}:{s}-{e}'.format(c=chr,s=start,e=end)
        out = '{o}/{s}.bcf'.format(o=tmpdir,s=start)
        #if __name__ == '__main__':
        f = pool.apply_async(mpileup_region, [region,out,bam_files])
        print (f)
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

def mpileup_gnuparallel(bam_files, ref, outpath, threads=4, callback=None):
    """Run mpileup in over multiple regions with GNU parallel, then concat vcf files.
    Assumes alignment to a bacterial reference with a single chromosome."""

    bam_files = ' '.join(bam_files)
    rawbcf = os.path.join(outpath,'raw.bcf')
    tmpdir = '/tmp'
    chr = tools.get_chrom(ref)
    length = tools.get_fasta_length(ref)
    x = np.linspace(1,length,threads+1,dtype=int)
    print (x)
    #split genome into blocks
    blocks=[]
    for i in range(len(x)):
        if i < len(x)-1:
            blocks.append((x[i],x[i+1]-1))

    outfiles = []
    regions = []
    for start,end in blocks:
        region = '"{c}":{s}-{e}'.format(c=chr,s=start,e=end)
        regions.append(region)
        out = '{o}/{s}.bcf'.format(o=tmpdir,s=start)
        outfiles.append(out)

    regstr = ' '.join(regions)
    filesstr = ' '.join(outfiles)
    annotatestr = '"AD,ADF,ADR,DP,SP,INFO/AD,INFO/ADF,INFO/ADR"'
    cmd = 'parallel bcftools mpileup -r {{1}} -a {a} -O b  --min-MQ 60 -o {{2}} -f {r} {b} ::: {reg} :::+ {o}'\
            .format(r=ref, reg=regstr, b=bam_files, o=filesstr, a=annotatestr)
    print (cmd)
    if callback != None:
        callback(cmd)
    subprocess.check_output(cmd, shell=True)
    #concat files
    cmd = 'bcftools concat {i} -O b -o {o}'.format(i=' '.join(outfiles),o=rawbcf)
    print (cmd)
    subprocess.check_output(cmd, shell=True)
    #remove temp files
    for f in outfiles:
        os.remove(f)
    return rawbcf

def variant_calling(bam_files, ref, outpath, relabel=True, threads=4,
                    callback=None, overwrite=False, filters=None, gff_file=None, **kwargs):
    """Call variants with bcftools"""

    st = time.time()
    if filters == None:
        filters = default_filter
    rawbcf = os.path.join(outpath,'raw.bcf')
    bcftoolscmd = tools.get_cmd('bcftools')
    if not os.path.exists(rawbcf) or overwrite == True:
        if platform.system() == 'Windows' or threads == 1:
            bam_files = ' '.join(bam_files)
            cmd = '{bc} mpileup -O b --min-MQ 60 -o {o} -f {r} {b}'.format(bc=bcftoolscmd,r=ref, b=bam_files, o=rawbcf)
            print (cmd)
            subprocess.check_output(cmd, shell=True)
        #if linux use mpileup in parallel to speed up
        else:
            #rawbcf = mpileup_multiprocess(bam_files, ref, outpath, threads=threads, callback=callback)
            rawbcf = mpileup_gnuparallel(bam_files, ref, outpath, threads=threads, callback=callback)

    #find snps
    vcfout = os.path.join(outpath,'calls.vcf')
    cmd = '{bc} call -V indels --ploidy 1 -m -v -o {v} {raw}'.format(bc=bcftoolscmd,v=vcfout,raw=rawbcf)
    if callback != None:
        callback(cmd)
    print (cmd)
    subprocess.check_output(cmd,shell=True)

    #relabel samples in vcf header
    if relabel == True:
        sample_file = os.path.join(outpath,'samples.txt')
        rlout = os.path.join(tempdir,'calls.vcf')
        cmd = 'bcftools reheader --samples {s} -o {o} {v}'.format(o=rlout,v=vcfout,s=sample_file)
        print(cmd)
        tmp = subprocess.check_output(cmd,shell=True)
        shutil.copy(rlout, vcfout)

    #filter variants
    final = os.path.join(outpath,'filtered.vcf.gz')
    cmd = '{bc} filter -i "{f}" -o {o} -O z {i}'.format(bc=bcftoolscmd,i=vcfout,o=final,f=filters)
    print (cmd)
    tmp = subprocess.check_output(cmd,shell=True)
    if callback != None:
        callback(cmd)

    #consequence calling
    if gff_file != None:
        csqout = os.path.join(outpath, 'csq.tsv')
        cmd = 'bcftools csq -f {r} -g {g} {f} -Ot -o {o}'.format(r=ref,g=gff_file,f=final,o=csqout)
        print (cmd)
        tmp = subprocess.check_output(cmd,shell=True)
        csqdf = read_csq_file(csqout)
        #get presence/absence matrix of csq mutations
        m = get_aa_snp_matrix(csqdf)
        m.to_csv(os.path.join(outpath,'csq.matrix'))
    print ('took %s seconds' %str(round(time.time()-st,0)))
    return final

def create_bam_labels(filenames):

    names = [os.path.basename(i).split('.')[0] for i in filenames]
    with open('samples.txt','w+') as file:
        for s in zip(bam_files,names):
            file.write('%s %s\n' %(s[0],s[1]))
    return

def trim_files(df, outpath, overwrite=False, threads=4, quality=30):
    """Batch trim fastq files"""

    method = 'cutadapt'
    if platform.system() == 'Windows':
        method = 'default'
    if not os.path.exists(outpath):
        os.makedirs(outpath, exist_ok=True)
    for i,row in df.iterrows():
        outfile = os.path.join(outpath, os.path.basename(row.filename))
        if not os.path.exists(outfile) or overwrite == True:
            tools.trim_reads(row.filename, outfile, threads=threads, quality=quality, method=method)
            print (outfile)
        df.loc[i,'trimmed'] = outfile
    return df

def read_csq_file(filename):
    """Read csq tsv outpt file into dataframe"""

    cols = ['1','sample','2','chrom','start','snp_type','gene','locus_tag','strand','feature_type','aa','nuc',]
    csqdf = pd.read_csv(filename,sep='[|\t]',comment='#',names=cols, engine='python')
    return csqdf

def get_aa_snp_matrix(df):
    """Get presence/absence matrix from csq calls table"""

    df=df.drop_duplicates(['gene','aa','sample'])
    x = df.set_index(['gene','aa','sample'])['start'].unstack('sample')
    x[x.notna()] = 1
    x = x.fillna(0)
    return x

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
        print ('options')
        print ('-------')
        for i in self.__dict__:
            print (i, ':', self.__dict__[i])
        print ()
        return

    def setup(self):
        """Setup main parameters"""

        if self.reference == None:
            self.reference = mbovis_genome
            self.gff_file = mbovis_gff
        self.filenames = get_files_from_paths(self.input)
        if self.threads == None:
            import multiprocessing
            self.threads = multiprocessing.cpu_count()
        else:
            self.threads = int(self.threads)
        df = get_samples(self.filenames, sep=self.labelsep)
        if len(df) == 0:
            print ('no samples provided')
            return
        df['read_length'] = df.filename.apply(tools.get_fastq_info)
        self.fastq_table = df
        sample_size = len(df['sample'].unique())
        print ('%s samples were loaded:' %sample_size)
        print ('----------------------')
        print (df)
        print ()
        s = check_samples_unique(df)
        if s == False:
            print ('samples names are not unique! try a different labelsep value.')
            return False
        print ('building index')
        aligners.build_bwa_index(self.reference)
        time.sleep(1)
        return True

    def run(self):
        """Run workflow"""

        #this master table tracks our outputs
        samples = self.fastq_table
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir, exist_ok=True)
        write_samples(samples, self.outdir)
        if len(samples)==0:
            print ('no samples found')
            return

        if self.trim == True:
            print ('trimming fastq files')
            print ('--------------------')
            trimmed_path = os.path.join(self.outdir, 'trimmed')
            samples = trim_files(samples, trimmed_path, self.overwrite,
                                  quality=self.quality, threads=self.threads)
            print ()
        print ('aligning files')
        print ('--------------')
        print ('Using reference genome: %s' %self.reference)
        path = os.path.join(self.outdir, 'mapped')
        samples = align_reads(samples, idx=self.reference, outdir=path,
                    threads=self.threads, overwrite=self.overwrite)
        print ()
        print ('calling variants')
        print ('----------------')
        bam_files = list(samples.bam_file.unique())
        self.vcf_file = variant_calling(bam_files, self.reference, self.outdir, threads=self.threads,
                                        gff_file=self.gff_file, filters=self.filters,
                                        overwrite=self.overwrite)
        print (self.vcf_file)
        print ()
        print ('making SNP matrix')
        print ('-----------------')
        snprecs, smat = tools.fasta_alignment_from_vcf(self.vcf_file)
        outfasta = os.path.join(self.outdir, 'core.fa')
        SeqIO.write(snprecs, outfasta, 'fasta')
        #write out sites matrix as txt file
        smat.to_csv(os.path.join(self.outdir,'core.txt'), sep=' ')
        print ()
        #save summary table
        summ = results_summary(samples)
        summ.to_csv(os.path.join(self.outdir,'summary.csv'),index=False)
        print ('Done. Sample summary:')
        print ('---------------------')
        print (summ)
        print ()

        if self.buildtree == True:
            print ('building tree')
            print ('-------------')
            if len(bam_files) <= 2:
                print ('Cannot build tree, too few samples.')
                return
            treefile = trees.run_RAXML(outfasta, outpath=self.outdir)
            if treefile == None:
                return
            print (treefile)
            #labelmap = dict(zip(sra.filename,sra.geo_loc_name_country))
            t,ts = trees.create_tree(treefile)#, labelmap)
            t.render(os.path.join(self.outdir, 'tree.png'))
        print ()
        return

def test_run():
    """Test run"""

    args = {'threads':8, 'outdir': 'testing', 'input':'mbovis_sra',
            'reference': None, 'overwrite':False}
    W = WorkFlow(**args)
    st = W.setup()
    W.run()
    return

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='snpgenie CLI tool. https://github.com/dmnfarrell/snpgenie')
    parser.add_argument("-i", "--input", action='append', dest="input", default=[],
                        help="input folder(s)", metavar="FILE")
    #parser.add_argument("-l", "--labels", dest="labels", default=[],
    #                    help="sample labels file, optional", metavar="FILE")
    parser.add_argument("-e", "--labelsep", dest="labelsep", default=',',
                        help="symbol to split the sample labels on")
    parser.add_argument("-r", "--reference", dest="reference", default=None,
                        help="reference genome filename", metavar="FILE")
    parser.add_argument("-g", "--gff", dest="gff_file", default=None,
                        help="reference gff, optional", metavar="FILE")
    parser.add_argument("-w", "--overwrite", dest="overwrite", action="store_true", default=False,
                        help="overwrite intermediate files")
    parser.add_argument("-m", "--trim", dest="trim", action="store_true", default=False,
                        help="whether to trim fastq files" )
    parser.add_argument("-q", "--quality", dest="quality", default=25,
                        help="trim quality" )
    parser.add_argument("-f", "--filters", dest="filters", default=default_filter,
                        help="variant calling post-filters" )
    parser.add_argument("-t", "--threads", dest="threads", default=None,
                        help="cpu threads to use")
    parser.add_argument("-b", "--buildtree", dest="buildtree", action="store_true", default=False,
                        help="whether to try to build a phylogenetic tree" )
    parser.add_argument("-o", "--outdir", dest="outdir",
                        help="Results folder", metavar="FILE")
    parser.add_argument("-v", "--version", dest="version", action="store_true",
                        help="Get version")
    parser.add_argument("-d", "--dummy", dest="dummy",  action="store_true",
                        default=False, help="Setup samples but don't run")

    args = vars(parser.parse_args())
    check_platform()
    print (datetime.datetime.now())
    #Log = Logger(os.path.join(args['outdir'], 'run.log'))
    #if args['test'] == True:
    #    test_run()
    if args['version'] == True:
        from . import __version__
        print ('snpgenie version %s' %__version__)
        print ('https://github.com/dmnfarrell/snpgenie')
    elif args['outdir'] == None:
        print ('No input or output folders provided. These are required.')
        print ('Example:')
        print ('snpgenie -r <reference> -i <input folder with fastq.gz files> -o <output folder>')
        print ('Use -h for more help on options.')
    else:
        W = WorkFlow(**args)
        st = W.setup()
        if st == True and args['dummy'] == False:
            W.run()

if __name__ == '__main__':
    main()
