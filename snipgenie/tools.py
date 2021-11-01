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

import sys,os,subprocess,glob,shutil,re,random,time
import platform
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Phylo, AlignIO
import numpy as np
import pandas as pd
#import pylab as plt
from gzip import open as gzopen

home = os.path.expanduser("~")
config_path = os.path.join(home,'.config','snipgenie')
bin_path = os.path.join(config_path, 'binaries')
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(module_path, 'data')

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """

    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)

def get_cmd(cmd):
    """Get windows version of a command if required"""

    if getattr(sys, 'frozen', False):
        cmd = resource_path('bin\%s.exe' %cmd)
    elif platform.system() == 'Windows':
        cmd = os.path.join(bin_path, '%s.exe' %cmd)
    return cmd

def move_files(files, path):
    if not os.path.exists(path):
        os.mkdir(path)
    for f in files:
        shutil.move(f, os.path.join(path,os.path.basename(f)))
    return

def get_attributes(obj):
    """Get non hidden and built-in type object attributes that can be persisted"""

    d={}
    import matplotlib
    allowed = [str,int,float,list,tuple,bool,matplotlib.figure.Figure]
    for key in obj.__dict__:
        if key.startswith('_'):
            continue
        item = obj.__dict__[key]
        if type(item) in allowed:
            d[key] = item
        elif type(item) is dict:
            if checkDict(item) == 1:
                d[key] = item
    return d

def set_attributes(obj, data):
    """Set attributes from a dict. Used for restoring settings in tables"""

    for key in data:
        try:
            obj.__dict__[key] = data[key]
        except Exception as e:
            print (e)
    return

def checkDict(d):
    """Check a dict recursively for non serializable types"""

    allowed = [str,int,float,list,tuple,bool]
    for k, v in d.items():
        if isinstance(v, dict):
            checkDict(v)
        else:
            if type(v) not in allowed:
                return 0
    return 1

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

def diffseqs(seq1,seq2):
    """Diff two sequences"""

    c=0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            c+=1
    return c

def snp_dist_matrix(aln):
    """Get pairwise snps distances from biopython
       Multiple Sequence Alignment object.
       returns: pandas dataframe
    """

    names=[s.id for s in aln]
    m = []
    for i in aln:
        x=[]
        for j in aln:
            d = diffseqs(i,j)
            x.append(d)
        #print (x)
        m.append(x)
    m = pd.DataFrame(m,index=names,columns=names)
    return m

def get_fasta_length(filename):
    """Get length of reference sequence"""

    from pyfaidx import Fasta
    refseq = Fasta(filename)
    key = list(refseq.keys())[0]
    l = len(refseq[key])
    return l

def get_chrom(filename):
    """Get chromosome name from fasta file"""

    rec = list(SeqIO.parse(filename, 'fasta'))[0]
    return rec.id

def get_fastq_info(filename):
    """Return fastq mean read length"""

    df = fastq_to_dataframe(filename)
    name = os.path.basename(filename).split('.')[0]
    rl = int(df.length.mean())
    return rl

def get_fastq_length(filename):
    """Return fastq number of reads"""

    cmd = 'expr $(zcat "%s" | wc -l) / 4' %filename
    tmp = subprocess.check_output(cmd, shell=True)
    l = int(tmp)
    return l

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

def make_blast_database(filename, dbtype='nucl'):
    """Create a blast db from fasta file"""

    cmd = get_cmd('makeblastdb')
    cline = '%s -dbtype %s -in %s' %(cmd,dbtype,filename)
    subprocess.check_output(cline, shell=True)
    return

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
                    compress=False, cmd='blastn', threads=4, show_cmd=False, **kwargs):
    """Blast a local database.
    Args:
        database: local blast db name
        query: sequences to query, list of strings or Bio.SeqRecords
    Returns:
        pandas dataframe with top blast results
    """

    if output == None:
        output = os.path.splitext(query)[0]+'_blast.txt'
    from Bio.Blast.Applications import NcbiblastnCommandline
    outfmt = '"6 qseqid sseqid qseq sseq pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore stitle"'
    cline = NcbiblastnCommandline(query=query, cmd=cmd, task='blastn', db=database,
                                 max_target_seqs=maxseqs,
                                 outfmt=outfmt, out=output,
                                 evalue=evalue, num_threads=threads, **kwargs)
    if show_cmd == True:
        print (cline)
    stdout, stderr = cline()
    return output

def blast_fasta(database, filename, **kwargs):
    """
    Blast a fasta file
    """

    remotedbs = ['nr','refseq_protein','pdb','swissprot']
    if database in remotedbs:
        remote_blast(database, filename, **kwargs)
    else:
        outfile = local_blast(database, filename, **kwargs)
    df = get_blast_results(filename=outfile)
    return df

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
    df = blast_file(database, 'tempseq.fa', **kwargs)
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
        rec = SeqRecord(Seq(row[seqkey]),id=str(row[idkey]),
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

def bam_to_fastq(filename):
    """bam to fastq using samtools"""

    samtoolscmd = get_cmd('samtools')
    name = os.path.basename(filename,threads=4)
    cmd = '{s} fastq -@ {t} {f} \
    -1 {n}_R1.fastq.gz -2 {n}_R2.fastq.gz \
    -0 /dev/null -s /dev/null -n'.format(s=samtoolscmd,f=filename,n=name,t=threads)
    print (cmd)
    subprocess.check_output(cmd, shell=True)
    return

def fastq_to_fasta(filename, out, size=1000):
    """Convert fastq to fasta
        size: limit to the first reads of total size
    """

    ext = os.path.splitext(filename)[1]
    if ext=='.gz':
        #fastq_parser = SeqIO.parse(gzopen(filename, "rt"), "fastq")
        #if gzip fails
        cmd = 'zcat "%s" | head -n %s > temp.fastq' %(filename, int(size))
        subprocess.check_output(cmd, shell=True)
        fastq_parser = SeqIO.parse(open('temp.fastq', "r"), "fastq")
    else:
        fastq_parser = SeqIO.parse(open(filename, "r"), "fastq")
    i=0
    recs=[]
    for fastq_rec in fastq_parser:
        i+=1
        if i>size:
            break
        recs.append(fastq_rec)
    SeqIO.write(recs, out, 'fasta')
    return

def records_to_dataframe(records, cds=False, nucl_seq=False):
    """Get features from a biopython seq record object into a dataframe
    Args:
        features: Bio SeqFeatures
        returns: a dataframe with a row for each cds/entry.
      """

    res = []
    for rec in records:
        featurekeys = []
        allfeat = []
        for (item, f) in enumerate(rec.features):
            x = f.__dict__
            quals = f.qualifiers
            x.update(quals)
            d = {}
            d['start'] = f.location.start
            d['end'] = f.location.end
            d['strand'] = f.location.strand
            d['id'] = rec.id
            d['feat_type'] = f.type
            for i in quals:
                if i in x:
                    if type(x[i]) is list:
                        d[i] = x[i][0]
                    else:
                        d[i] = x[i]
            allfeat.append(d)
            if nucl_seq == True:
                nseq = str(rec.seq)[f.location.start:f.location.end]
                if f.location.strand == -1:
                    d['sequence'] = str(Seq(nseq).reverse_complement())
                else:
                    d['sequence'] = nseq
        quals = list(quals.keys())+['id','start','end','strand','feat_type','sequence']
        df = pd.DataFrame(allfeat,columns=quals)
        if 'translation' in df.keys():
            df['length'] = df.translation.astype('str').str.len()
        res.append(df)

    if len(df) == 0:
        print ('ERROR: genbank file return empty data, check that the file contains protein sequences '\
               'in the translation qualifier of each protein feature.' )
    return pd.concat(res)

def genbank_to_dataframe(infile, cds=False):
    """Get genome records from a genbank file into a dataframe
      returns a dataframe with a row for each cds/entry"""

    recs = list(SeqIO.parse(infile,'genbank'))
    df = records_to_dataframe(recs, cds, nucl_seq=True)
    return df

def gff_to_records(gff_file):
    """Get features from gff file"""

    if gff_file is None or not os.path.exists(gff_file):
        return
    from BCBio import GFF
    in_handle = open(gff_file,'r')
    recs = list(GFF.parse(in_handle))
    in_handle.close()
    return recs

def features_summary(df):
    """SeqFeatures dataframe summary"""

    def hypo(val):
        val = val.lower()
        kwds=['hypothetical','conserved protein','unknown protein']
        for k in kwds:
            if k in val:
                return True
        else:
            return False
    coding = df[df.feat_type=='CDS']
    trna = df[df.feat_type=='tRNA']
    products = coding[coding['product'].notnull()]
    cdstrans = coding[coding.translation.notnull()]
    hypo = products[products['product'].apply(hypo)]
    #pseudo = df[ (df.feat_type == 'gene') & (df.pseudo.notnull())]
    notags = df[df.locus_tag.isnull()]
    repeats = df[ (df.feat_type == 'repeat_region')]
    s = {}
    s['total features'] = len(df)
    s['coding sequences'] = len(coding)
    s['cds with translations'] = len(cdstrans)
    s['cds with gene names'] = len(coding.gene.dropna())
    s['hypothetical'] = len(hypo)
    #s['pseudogenes'] = len(pseudo)
    s['trna'] = len(trna)
    s['repeat_region'] = len(repeats)
    s['no locus tags'] =  len(notags)
    if len(cdstrans)>0:
        avlength = int(np.mean([len(i) for i in cdstrans.translation]))
        s['mean sequence length'] =  avlength

    return s

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

def vcf_to_dataframe(vcf_file):
    """
    Convert a multi sample vcf to dataframe. Records each samples FORMAT fields.
    Args:
        vcf_file: input multi sample vcf
    Returns: pandas DataFrame
    """

    import vcf
    ext = os.path.splitext(vcf_file)[1]
    if ext == '.gz':
        file = gzopen(vcf_file, "rt")
    else:
        file = open(vcf_file)
    vcf_reader = vcf.Reader(file,'r')
    res=[]
    cols = ['sample','REF','ALT','mut','DP','ADF','ADR','AD','chrom','var_type',
            'sub_type','start','end','QUAL']
    i=0
    for rec in vcf_reader:
        #if i>50:
        #    break
        x = [rec.CHROM, rec.var_type, rec.var_subtype, rec.start, rec.end, rec.QUAL]
        for sample in rec.samples:
            if sample.gt_bases == None:
                mut=''
                row = [sample.sample, rec.REF, sample.gt_bases, mut, 0,0,0,0]
            elif rec.REF != sample.gt_bases:
                mut = str(rec.end)+rec.REF+'>'+sample.gt_bases
                cdata = sample.data
                row = [sample.sample, rec.REF, sample.gt_bases, mut, cdata[2], cdata[4] ,cdata[5], cdata[6]] + x
            else:
                mut = str(rec.end)+rec.REF
                #inf = sample.site.INFO
                cdata = sample.data
                row = [sample.sample, rec.REF, sample.gt_bases, mut, cdata[2], cdata[4] ,cdata[5], cdata[6]] + x

            res.append(row)
    res = pd.DataFrame(res,columns=cols)
    res = res[~res.start.isnull()]
    #res['start'] = res.start.astype(int)
    #res['end'] = res.end.astype(int)
    return res

def get_snp_matrix(df):
    """SNP matrix from multi sample vcf dataframe"""

    df = df.drop_duplicates(['mut','sample'])
    x = df.set_index(['mut','sample'])['start'].unstack('sample')
    x[x.notna()] = 1
    x = x.fillna(0)
    return x

def plot_fastq_qualities(filename, ax=None, limit=10000):
    """Plot fastq qualities"""

    if not os.path.exists(filename):
        return
    import matplotlib.patches as patches
    import pylab as plt
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
    n = int(len(df.columns)/50)
    #df = df[df.columns[::n]]
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
    n = int(l/10)
    ax.set_xticks(np.arange(0, l, n))
    ax.set_xticklabels(np.arange(0, l, n))

    ax.set_xlabel('position(bp)')
    ax.set_xlim((0,l))
    ax.set_ylim((0,40))
    ax.set_title('per base sequence quality')
    return

def normpdf(x, mean, sd):
    """Normal distribution function"""

    import math
    var = float(sd)**2
    denom = (2*math.pi*var)**.5
    num = math.exp(-(float(x)-float(mean))**2/(2*var))
    return num/denom

def plot_fastq_gc_content(filename, ax=None, limit=50000):
    """Plot fastq gc conent"""

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

def fastq_quality_report(filename, figsize=(7,5), **kwargs):
    """Fastq quality plots"""

    import pylab as plt
    fig,ax = plt.subplots(2,1, figsize=figsize, dpi=100, facecolor=(1,1,1), edgecolor=(0,0,0))
    axs = ax.flat
    plot_fastq_qualities(filename, ax=axs[0], **kwargs)
    plot_fastq_gc_content(filename, ax=axs[1])
    plt.tight_layout()
    return fig

def pdf_qc_reports(filenames, outfile='qc_report.pdf'):
    """Save pdf reports of fastq file quality info"""

    import pylab as plt
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(outfile) as pdf:
        for f in sorted(filenames):
            print (f)
            fig = fastq_quality_report(f)
            fig.subplots_adjust(top=.9)
            fig.suptitle(os.path.basename(f))
            pdf.savefig(fig)
            plt.clf()

def concat_seqrecords(recs):
    """Join seqrecords together"""
    concated = Seq("")
    for r in recs:
        concated += r.seq
    return SeqRecord(concated, id=recs[0].id)

def fasta_alignment_from_vcf(vcf_file, callback=None, uninformative=False, omit=None):
    """
    Get core SNP site calls as sequences from a multi sample vcf file.
    Args:
        vcf_file: multi-sample vcf (e.g. produced by app.variant_calling)
        uninformative: whether to include uninformative sites
        omit: list of samples to exclude if required
    """

    import vcf
    from collections import defaultdict
    vcf_reader = vcf.Reader(open(vcf_file, 'rb'))
    #print (vcf_reader.samples)
    def default():
        return []
    result = defaultdict(default)
    sites = []
    result['ref'] = []
    missing = []
    unf = []
    for record in vcf_reader:
        S = {sample.sample: sample.gt_bases for sample in record.samples}
        if omit != None:
            for o in omit:
                del S[o]
        #if any missing samples at the site we don't add
        if None in S.values():
            #for i in S:
            #    if S[i]==None:
            #        print (i)
            missing.append(record.POS)
            continue
        if uninformative == False:
            u = set(S.values())
            if len(u) == 1:
                unf.append(record.POS)
                continue
        result['ref'].append(record.REF)
        for name in S:
            result[name].append(S[name])
        sites.append(record.POS)

    print ('found %s sites' %len(sites))
    print ('%s sites with at least one missing sample' %len(missing))
    if uninformative == False:
        print ('%s uninformative sites' %len(unf))
    if len(sites)==0:
        print ('no sites found may mean one sample is too different')
    recs = []
    for sample in result:
        seq = ''.join(result[sample])
        seqrec = SeqRecord(Seq(seq),id=sample)
        recs.append(seqrec)
        #print (len(seqrec))

    smat = pd.DataFrame(result)
    smat.index = sites
    smat.index.rename('pos', inplace=True)
    return recs, smat

def samtools_flagstat(filename):
    """Parse samtools flagstat output into dictionary"""

    samtoolscmd = get_cmd('samtools')
    cmd = '{s} flagstat {f}'.format(f=filename,s=samtoolscmd)
    tmp = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    x = tmp.split('\n')
    x = [int(i.split('+')[0]) for i in x[:-1]]
    #print (x)
    cols = ['total','secondary','supplementary','duplicates','mapped',
            'paired','read1','read2','properly paired','with itself','singletons']
    d = {}
    for c,v in zip(cols,x):
        d[c] = v
    return d

def samtools_tview(bam_file, chrom, pos, width=200, ref='', display='T'):
    """View bam alignment with samtools"""

    samtoolscmd = get_cmd('samtools')
    os.environ["COLUMNS"] = str(width)
    cmd = '{sc} tview {b} -p {c}:{p} -d {d} {r}'\
            .format(b=bam_file,c=chrom,p=pos,d=display,r=ref,sc=samtoolscmd)
    #print (cmd)
    tmp = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    return tmp

def samtools_depth(bam_file, chrom=None, start=None, end=None):
    """Get depth from bam file"""

    samtoolscmd = get_cmd('samtools')
    if chrom != None and start != None:
        cmd = '{sc} depth -r {c}:{s}-{e} {b}'.format(b=bam_file,c=chrom,s=start,e=end,sc=samtoolscmd)
    else:
        cmd = '{sc} depth {b}'.format(b=bam_file,sc=samtoolscmd)
    tmp=subprocess.check_output(cmd, shell=True, universal_newlines=True)
    from io import StringIO
    c = pd.read_csv(StringIO(tmp),sep='\t',names=['chr','pos','depth'])
    return c

def get_mean_depth(bam_file, chrom=None, start=None, end=None, how='mean'):
    """Get mean depth from bam file"""

    c = samtools_depth(bam_file, chrom, start, end)
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

def gff_bcftools_format(in_file, out_file):
    """Convert a genbank file to a GFF format that can be used in bcftools csq.
    see https://github.com/samtools/bcftools/blob/develop/doc/bcftools.txt#L1066-L1098.
    Args:
        in_file: genbank file
        out_file: name of GFF file
    """

    from BCBio import GFF
    in_handle = open(in_file)
    out_handle = open(out_file, "w")
    from Bio.SeqFeature import SeqFeature
    from Bio.SeqFeature import FeatureLocation
    from copy import copy

    #recs = GFF.parse(in_handle)
    recs = SeqIO.parse(in_file,format='gb')
    for record in recs:
        #make a copy of the record as we will be changing it during the loop
        new = copy(record)
        new.features = []
        #loop over all features
        for i in range(0,len(record.features)):
            feat = record.features[i]
            q = feat.qualifiers
            if not 'locus_tag' in q:
                continue
            #print (q)
            #remove some unecessary qualifiers
            for label in ['note','translation','product','experiment']:
                if label in q:
                    del q[label]
            if(feat.type == "CDS"):
                #use the CDS feature to create the new lines
                tag = q['locus_tag'][0]
                q['ID'] = 'CDS:%s' %tag
                q['Parent'] = 'transcript:%s' %tag
                q['biotype'] = 'protein_coding'
                #create mRNA feature
                m = SeqFeature(feat.location,type='mRNA',strand=feat.strand)
                q = m.qualifiers
                q['ID'] = 'transcript:%s' %tag
                q['Parent'] = 'gene:%s' %tag
                q['biotype'] = 'protein_coding'
                new.features.append(m)
            elif(record.features[i].type == "gene"):
                #edit the gene feature
                q=feat.qualifiers
                q['ID'] = 'gene:%s' %q['locus_tag'][0]
                q['biotype'] = 'protein_coding'
                if 'gene' in q:
                    q['Name'] = q['gene']
            new.features.append(feat)
        #write the new features to a GFF
        GFF.write([new], out_handle)
        return

def get_spoligotype(filename, reads_limit=500000, threshold=2):
    """Get mtb spoligotype from WGS reads"""

    ref = os.path.join(datadir, 'dr_spacers.fa')
    #convert reads to fasta
    fastq_to_fasta(filename, 'temp.fa', reads_limit)
    #make blast db from reads
    make_blast_database('temp.fa')
    #blast spacers to db
    bl = blast_fasta('temp.fa', ref, evalue=0.1,
                           maxseqs=100000, show_cmd=False)
    bl=bl[(bl.qcovs>95) & (bl.mismatch<2)]
    x = bl.groupby('qseqid').agg({'pident':np.size}).reset_index()
    x = x[x.pident>=threshold]
    found = list(x.qseqid)
    s=[]
    for i in range(1,44):
        if i in found:
            s.append('1')
        else:
            s.append('0')
    s =''.join(s)
    #print (s)
    return s

def get_sb_number(binary_str):
    """Get SB number from binary pattern usinf database reference"""
    
    df = pd.read_csv(os.path.join(datadir, 'Mbovis.org_db.csv'))
    x = df[df['binary'] == binary_str]
    if len(x) == 0:
        return
    else:
        return x.iloc[0].SB
