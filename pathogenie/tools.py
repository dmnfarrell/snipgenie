"""
    Various methods for bacterial genomics.
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

home = os.path.expanduser("~")
config_path = os.path.join(home,'.config','pathogenie')
bin_path = os.path.join(config_path, 'binaries')

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """

    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)

def get_cmd(cmd):
    """Get windows version of a command if required"""

    if getattr(sys, 'frozen', False):
        cmd = tools.resource_path('bin/%s.exe' %cmd)
    elif platform.system() == 'Windows':
        cmd = os.path.join(bin_path, '%s.exe' %cmd)
    return cmd

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.__next__()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

def get_fasta_length(filename):
    """Get length of reference sequence"""

    from pyfaidx import Fasta
    refseq = Fasta(filename)
    key = list(refseq.keys())[0]
    l = len(refseq[key])
    return l

def get_chrom(filename):
    rec = list(SeqIO.parse(filename, 'fasta'))[0]
    return rec.id

def get_fastq_info(filename):

    df = fastq_to_dataframe(filename)
    name = os.path.basename(filename).split('.')[0]
    rl = int(df.length.mean())
    return rl

def align_info(bamfile):

    cmd = 'samtools flagstat %s' %bamfile
    temp=subprocess.check_output(cmd, shell=True)
    print (temp)
    return

def clustal_alignment(filename=None, seqs=None, command="clustalw"):
    """Align 2 sequences with clustal"""

    if filename == None:
        filename = 'temp.faa'
        SeqIO.write(seqs, filename, "fasta")
    name = os.path.splitext(filename)[0]
    from Bio.Align.Applications import ClustalwCommandline
    cline = ClustalwCommandline(command, infile=filename)
    stdout, stderr = cline()
    align = AlignIO.read(name+'.aln', 'clustal')
    return align

def get_blast_results(filename):
    """
    Get blast results into dataframe. Assumes column names from local_blast method.
    Returns:
        dataframe
    """

    cols = ['qseqid','sseqid','qseq','sseq','pident','qcovs','length','mismatch','gapopen',
            'qstart','qend','sstart','send','evalue','bitscore','stitle']
    res = pd.read_csv(filename, names=cols, sep='\t')
    #res = res[res['pident']>=ident]
    return res

def local_blast(database, query, output=None, maxseqs=50, evalue=0.001,
                    compress=False, cmd='blastn', cpus=4, show_cmd=False, **kwargs):
    """Blast a local database.
    Args:
        database: local blast db name
        query: sequences to query, list of strings or Bio.SeqRecords
    Returns:
        pandas dataframe with top blast results
    """

    if output == None:
        output = os.path.splitext(query)[0]+'_blast.txt'
    from Bio.Blast.Applications import NcbiblastxCommandline
    outfmt = '"6 qseqid sseqid qseq sseq pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore stitle"'
    cline = NcbiblastxCommandline(query=query, cmd=cmd, db=database,
                                 max_target_seqs=maxseqs,
                                 outfmt=outfmt, out=output,
                                 evalue=evalue, num_threads=cpus, **kwargs)
    if show_cmd == True:
        print (cline)
    stdout, stderr = cline()
    return

def blast_sequences(database, seqs, labels=None, **kwargs):
    """
    Blast a set of sequences to a local or remote blast database
    Args:
        database: local or remote blast db name
                  'nr', 'refseq_protein', 'pdb', 'swissprot' are valide remote dbs
        seqs: sequences to query, list of strings or Bio.SeqRecords
        labels: list of id names for sequences, optional but recommended
    Returns:
        pandas dataframe with top blast results
    """

    remotedbs = ['nr','refseq_protein','pdb','swissprot']
    res = []
    if not type(seqs) is list:
        seqs = [seqs]
    if labels is None:
        labels = seqs
    recs=[]
    #print (labels)
    for seq, name in zip(seqs,labels):
        if type(seq) is not SeqRecord:
            rec = SeqRecord(Seq(seq),id=name)
        else:
            rec = seq
            name = seq.id
        recs.append(rec)
    SeqIO.write(recs, 'tempseq.fa', "fasta")
    if database in remotedbs:
        remote_blast(database, 'tempseq.fa', **kwargs)
    else:
        local_blast(database, 'tempseq.fa', **kwargs)
    df = get_blast_results(filename='tempseq_blast.txt')
    return df

def dataframe_to_fasta(df, seqkey='translation', idkey='locus_tag',
                     descrkey='description',
                     outfile='out.faa'):
    """Genbank features to fasta file"""

    seqs=[]
    for i,row in df.iterrows():
        if descrkey in df.columns:
            d=row[descrkey]
        else:
            d=''
        rec = SeqRecord(Seq(row[seqkey]),id=row[idkey],
                            description=d)
        seqs.append(rec)
    SeqIO.write(seqs, outfile, "fasta")
    return outfile

def fasta_to_dataframe(infile, header_sep=None, key='name', seqkey='sequence'):
    """Get fasta proteins into dataframe"""

    recs = SeqIO.parse(infile,'fasta')
    keys = [key,seqkey,'description']
    data = [(r.name,str(r.seq),str(r.description)) for r in recs]
    df = pd.DataFrame(data,columns=(keys))
    df['type'] = 'CDS'
    #fix bad names
    if header_sep not in ['',None]:
        df[key] = df[key].apply(lambda x: x.split(header_sep)[0],1)
    df[key] = df[key].str.replace('|','_')
    return df

def fastq_to_dataframe(filename, size=1000):
    """Convert fastq to dataframe.
        size: limit to the first reads of total size
        Returns: dataframe with reads
    """

    ext = os.path.splitext(filename)[1]
    if ext=='.gz':
        fastq_parser = SeqIO.parse(gzopen(filename, "rt"), "fastq")
    else:
        fastq_parser = SeqIO.parse(open(filename, "r"), "fastq")
    i=0
    res=[]
    for fastq_rec in fastq_parser:
        #print (fastq_rec.seq)
        i+=1
        if i>size:
            break
        res.append([fastq_rec.id, str(fastq_rec.seq)])
    df = pd.DataFrame(res, columns=['id','seq'])
    df['length'] = df.seq.str.len()
    return df

def features_to_dataframe(features, cds=False, id=''):
    """Get features from a biopython seq record object into a dataframe
    Args:
        features: bio seqfeatures
       returns: a dataframe with a row for each cds/entry.
      """

    featurekeys = []
    allfeat = []
    for (item, f) in enumerate(features):
        x = f.__dict__
        quals = f.qualifiers
        x.update(quals)
        d = {}
        d['start'] = f.location.start
        d['end'] = f.location.end
        d['strand'] = f.location.strand
        d['id'] = id
        d['feat_type'] = f.type
        for i in quals:
            if i in x:
                if type(x[i]) is list:
                    d[i] = x[i][0]
                else:
                    d[i] = x[i]
        allfeat.append(d)
    quals = list(quals.keys())+['id','start','end','strand','feat_type']
    df = pd.DataFrame(allfeat,columns=quals)
    if 'translation' in df.keys():
        df['length'] = df.translation.astype('str').str.len()

    if len(df) == 0:
        print ('ERROR: genbank file return empty data, check that the file contains protein sequences '\
               'in the translation qualifier of each protein feature.' )
    return df

def gff_to_features(gff_file):
    """Get features from gff file"""

    if gff_file is None or not os.path.exists(gff_file):
        return
    from BCBio import GFF
    in_handle = open(gff_file,'r')
    rec = list(GFF.parse(in_handle))[0]
    in_handle.close()
    return rec.features

def trim_reads_default(filename,  outfile, right_quality=35):
    """Trim adapters - built in method"""

    #trimmed = []
    fastq_parser = SeqIO.parse(gzopen(filename, "rt"), "fastq")
    c=0
    out = gzopen(outfile, "wt")
    for record in fastq_parser:
        score = record.letter_annotations["phred_quality"]
        for i in range(len(score)-1,0,-1):
            if score[i] >= right_quality:
                break
        #trimmed.append(record[:i])

        SeqIO.write(record[:i],out,'fastq')
    return

def trim_reads(filename, outfile, adapter=None, quality=20,
                method='cutadapt', threads=4):
    """Trim adapters using cutadapt"""

    #if adapter is not None and not type(adapter) is str:
    #    print ('not valid adapter')
    #    return
    if method == 'default':
        trim_reads_default(filename,  outfile, right_quality=quality)
    elif method == 'cutadapt':
        if adapter != None:
            cmd = 'cutadapt -O 5 -q {q} -a {a} -j {t} {i} -o {o}'.format(a=adapter,i=filename,o=outfile,t=threads,q=quality)
        else:
            cmd = 'cutadapt -O 5 -q {q} {i} -j {t} -o {o}'.format(i=filename,o=outfile,t=threads,q=quality)
        print (cmd)
        result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    return

def vcf_to_dataframe(vcf_file, quality=30):

    import vcf
    ext = os.path.splitext(vcf_file)[1]
    if ext == '.gz':
        file = gzopen(vcf_file, "rt")
    else:
        file = open(vcf_file)
    vcf_reader = vcf.Reader(file,'r')
    #print (vcf_reader.filters)
    res=[]
    cols = ['chrom','var_type','sub_type','start','end','REF','ALT','QUAL','DP']
    for rec in vcf_reader:
        x = [rec.CHROM, rec.var_type, rec.var_subtype, rec.start, rec.end, rec.REF, str(rec.ALT[0]),
            rec.QUAL, rec.INFO['DP']]
        #print (rec.__dict__)
        #print (rec.INFO.keys())
        #for call in rec.samples:
        #    print (call.sample, call.data, rec.genotype(call.sample))
        res.append(x)
        #print (x)
    res = pd.DataFrame(res,columns=cols)
    return res

def plot_fastq_qualities(filename, ax=None, limit=10000):

    import matplotlib.patches as patches
    fastq_parser = SeqIO.parse(gzopen(filename, "rt"), "fastq")
    res=[]
    c=0
    for record in fastq_parser:
        score=record.letter_annotations["phred_quality"]
        res.append(score)
        c+=1
        if c>limit:
            break
    df = pd.DataFrame(res)
    l = len(df.T)+1

    if ax==None:
        f,ax=plt.subplots(figsize=(12,5))
    rect = patches.Rectangle((0,0),l,20,linewidth=0,facecolor='r',alpha=.4)
    ax.add_patch(rect)
    rect = patches.Rectangle((0,20),l,8,linewidth=0,facecolor='yellow',alpha=.4)
    ax.add_patch(rect)
    rect = patches.Rectangle((0,28),l,12,linewidth=0,facecolor='g',alpha=.4)
    ax.add_patch(rect)
    df.mean().plot(ax=ax,c='black')
    boxprops = dict(linestyle='-', linewidth=1, color='black')
    df.plot(kind='box', ax=ax, grid=False, showfliers=False,
            color=dict(boxes='black',whiskers='black')  )
    ax.set_xticks(np.arange(0, l, 5))
    ax.set_xticklabels(np.arange(0, l, 5))
    ax.set_xlabel('position(bp)')
    ax.set_xlim((0,l))
    ax.set_ylim((0,40))
    ax.set_title('per base sequence quality')
    return

def normpdf(x, mean, sd):
    import math
    var = float(sd)**2
    denom = (2*math.pi*var)**.5
    num = math.exp(-(float(x)-float(mean))**2/(2*var))
    return num/denom

def plot_fastq_gc_content(filename, ax=None, limit=50000):

    from Bio.SeqUtils import GC
    if ax==None:
        f,ax=plt.subplots(figsize=(12,5))
    df = fastq_to_dataframe(filename, size=limit)
    gc = df.seq.apply(lambda x: GC(x))
    gc.hist(ax=ax,bins=150,color='black',grid=False,histtype='step',lw=2)
    ax.set_xlim((0,100))
    x=np.arange(1,100,.1)
    f = [normpdf(i, gc.mean(), gc.std()) for i in x]
    ax2=ax.twinx()
    ax2.plot(x,f)
    ax2.set_ylim(0,max(f))
    ax.set_title('GC content',size=15)
    return

def fasta_alignment_from_vcf(vcf_file, ref, callback=None):
    """Get a fasta alignment for all snp sites in a multi sample
    vcf file.
    Args:
        vcf_file: input vcf
        ref: the reference sequence
        callback: optional function to direct output
    """

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

    sites_matrix = {}
    sites_matrix['ref'] = list(refrec)
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
        sites_matrix[sample] = list(seqrec)

    #smat is a dataframe matrix of the positions and genotype
    smat = pd.DataFrame(sites_matrix)
    smat.index = sites
    return result, smat

def get_bam_depth(filename, how='mean'):
    """Get depth from bam file"""

    cmd = 'samtools depth %s ' %filename
    tmp=subprocess.check_output(cmd, shell=True)
    from io import StringIO
    c = pd.read_csv(StringIO(tmp.decode()),sep='\t',names=['chr','pos','depth'])
    if how == 'mean':
        return c.depth.mean().round(2)
    else:
        return c.depth.sum().round(2)

def fetch_sra_reads(df,path):
    """Download a set of reads from SRA using dataframe with runs"""

    for i,r in df.iterrows():
        files = glob.glob(os.path.join(path,r.Run+'*'))
        if len(files)==0:
            cmd = 'fastq-dump --split-files {n} --outdir {o}'.format(n=r.Run,o=path)
            print (cmd)
            subprocess.check_output(cmd,shell=True)
