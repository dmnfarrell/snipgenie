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
import hashlib
from io import StringIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Phylo, AlignIO, Align
import numpy as np
import pandas as pd
import pylab as plt
import matplotlib as mpl
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
            if check_dict(item) == 1:
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

def gunzip(infile, outfile):
    """Gunzip a file"""

    import gzip
    import shutil
    with gzip.open(infile, 'rb') as f_in:
        with open(outfile, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return

def check_dict(d):
    """Check a dict recursively for non serializable types"""

    allowed = [str,int,float,list,tuple,bool]
    for k, v in d.items():
        if isinstance(v, dict):
            check_dict(v)
        else:
            if type(v) not in allowed:
                return 0
    return 1

def is_folder_empty(folder_path):
    """
    Checks if a folder is empty or not.

    Parameters:
        folder_path (str): The path to the folder.

    Returns:
        bool: True if the folder is empty, False otherwise.
    """
    try:
        # Check if the path exists and is a directory
        if not os.path.isdir(folder_path):
            raise ValueError(f"The path '{folder_path}' is not a valid directory.")

        # Check if the folder is empty
        return len(os.listdir(folder_path)) == 0
    except Exception as e:
        print(f"Error: {e}")
        return False

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

def random_hex_color():
    """random hex color"""

    r = lambda: np.random.randint(0,255)
    c='#%02X%02X%02X' % (r(),r(),r())
    return c

def random_hex_colors(n=1,seed=None):
    if seed != None:
        np.random.seed(seed)
    return [random_hex_color() for i in range(n)]

def colormap_colors(colormap_name, n):
    """Colors list from mpl colormap"""

    colormap = plt.cm.get_cmap(colormap_name, n)
    colors = [mpl.colors.rgb2hex(colormap(i)) for i in range(n)]
    return colors

def colormap_from_labels(colormap_name, labels):
    """Get dict of colors mapping to labels using mpl colormap"""

    n = len(labels)
    colormap = plt.cm.get_cmap(colormap_name, n)
    colors = {labels[i]: mpl.colors.rgb2hex(colormap(i)) for i in range(n)}
    return colors

def df_html(df, fontsize='8pt'):
    """Create df html enclosed in some css"""

    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>Styled DataFrame</title>
        <style>
            body {{
                font-family: monospace;
                font-size: {fontsize};
                margin: 0;
                padding: 0;
            }}
            table {{
                border-collapse: collapse;
                width: 100%;
            }}
            th, td {{
                border: 1px solid #ddd;
                padding: 8px;
            }}
            th {{
                background-color: white;
                text-align: left;
                position: sticky;
                top: 0;
            }}

        </style>
    </head>
    <body>
        {df.to_html()}
    </body>
    </html>
    """
    return html

def diffseqs(seq1,seq2):
    """Diff two sequences"""

    c=0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            c+=1
    return c

def compute_snp_count(args):
    i, j, seq1, seq2 = args
    return i, j, np.sum(np.fromiter((c1 != c2 for c1, c2 in zip(seq1, seq2)), dtype=int))

def snp_dist_matrix(aln, threads=4):
    """
    Compute the number of Single Nucleotide Polymorphisms (SNPs)
    between sequences in a Biopython alignment.
    Args:
        aln:
            Biopython multiple sequence alignment object.
    Returns:
        A matrix as pandas dataframe.
    """

    from multiprocessing import Pool
    names = [s.id for s in aln]
    num_sequences = len(aln)
    matrix = np.zeros((num_sequences, num_sequences), dtype=int)
    sequences = [str(s.seq) for s in aln]

    args_list = []
    for i in range(num_sequences):
        seq1 = sequences[i]
        for j in range(i + 1, num_sequences):
            seq2 = sequences[j]
            args_list.append((i, j, seq1, seq2))

    with Pool(processes=threads) as pool:
        results = pool.map(compute_snp_count, args_list)

    for i, j, snp_count in results:
        matrix[i, j] = snp_count
        matrix[j, i] = snp_count

    m = pd.DataFrame(matrix, index=names, columns=names)
    return m

def update_snp_dist_matrix(aln, snpdist=None):
    """
    Compute a new SNP distance matrix given an updated
    alignment and previous matrix. The alignment should contain
    ALL the samples in the previous matrix. Designed so we can add
    one or few samples without recomputing all values.
    Args:
        aln: Biopython multiple sequence alignment object
        snpdist: current snp distance matrix to be updated
    returns:
        a matrix as pandas dataframe
    """

    present = snpdist.index
    #print (present)
    new_aln = [s for s in aln if s.id not in present]
    #print (len(new_aln))
    new=[s.id for s in new_aln]
    print ('%s samples to be added' %len(new))
    #create the new matrix
    matrix = snpdist.copy()
    new_index = snpdist.columns.tolist() + new
    matrix = snpdist.reindex(new_index)
    matrix = matrix.reindex(new_index, axis=1)

    for s1 in aln:
        seq1 = str(s1.seq)
        x=s1.id
        #print (x)
        for s2 in aln:
            y=s2.id
            if not np.isnan(matrix.loc[x,y]):
                continue
            else:
                seq2 = str(s2.seq)
                snp_count = sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
                matrix.loc[x,y] = snp_count
                matrix.loc[x,y] = snp_count

    #print (matrix[matrix.ref.isnull()])
    matrix = matrix.astype(int)
    return matrix

def alignment_from_snps(df):
    """MultipleSeqAlignment object from core snps"""

    if 'pos' in df.columns:
        df = df.set_index('pos')
    #transpose
    df = df.T
    seqs=[]
    for i,r in df.iterrows():
        s = ''.join(list(r))
        seqs.append(SeqRecord(Seq(s),id=i))

    aln = Align.MultipleSeqAlignment(seqs)
    return aln

def combine_core_snps(core, new):
    """
    Add new sample(s) to a core snps file by combining both tables.
    Args:
        core: dataframe of core snps from core.txt
        new: dataframe of new snps from core.txt
    returns:
        dataframe with combined snps
    """

    #new=new.drop(columns='ref')
    ncols = new.columns[1:]
    print (ncols)
    df = core.merge(new,on='pos',how='outer')
    #print (df)
    #fill first empty cols with ref
    for col in ncols:
        df[col] = df[col].fillna(df['ref_x'])
    currcols = core.columns
    #fill remainder
    df = df.apply(lambda x: x.fillna(df.ref_y))
    print ('%s->%s SNPs added' %(len(core),len(df)))
    df['ref'] = df.ref_y.fillna(df.ref_x)
    #print (df[df.ref.isnull()])
    df=df.drop(columns=['ref_y','ref_x'])
    return df

def get_closest_samples(distance_matrix, row_index, n):
    """
    Get the closest N samples for a given row in a distance matrix DataFrame.
    """
    row = distance_matrix.loc[row_index]
    closest = row.nsmallest(n + 1)[1:]
    return closest

def get_within_distance(distance_matrix, row_index, n):
    row = distance_matrix.loc[row_index]
    x = row[row<n]
    return x

def dist_matrix_to_mst(distance_matrix, df=None, colorcol=None, labelcol=None, colormap=None,
                       cmap_name='Set1', node_size=4, font_size=6, with_labels=False,
                       edge_labels=False, legend_loc=(1, .7), ax=None):
    """
    Dist matrix to minimum spanning tree
    Args:
        distance_matrix: matrix as dataframe
        df: meta data with index corresponding to node names of tree
        colorcol: column in meta table to color nodes by
        labelcol: column in meta table to label nodes by
        colormap: a mapping of node names to colors
        cmap_name: name of an mpl colormap to apply instead of providing colormap
        node_size: size of nodes
        with_labels: whether to plot labels on nodes
        legend_loc: location of legend
    """

    from . import plotting
    import pylab as plt
    if ax == None:
        fig,ax=plt.subplots()
    import networkx as nx
    from networkx.drawing.nx_agraph import graphviz_layout

    G = nx.Graph()

    for i, row in distance_matrix.iterrows():
        for j, weight in row.items():
            G.add_edge(i, j, weight=weight)

    T = nx.minimum_spanning_tree(G, algorithm='kruskal')
    # Compute edge lengths based on distances
    edge_lengths = [T[u][v]['weight'] for u, v in T.edges()]
    # Plot the minimum spanning tree with edge lengths proportional to distances
    pos = graphviz_layout(T)
    labels = nx.get_edge_attributes(T, 'weight')
    if df is not None:
        l = [label for label in T.nodes if label in df.index]
        df = df.loc[l]
        if colormap is None:
            colors,cmap = plotting.get_color_mapping(df, colorcol, cmap_name)
        else:
            #custom colormap if provided
            colors = [colormap[i] if i in colormap else 'black' for i in df[colorcol]]
            cmap = colormap
        #print (cmap)
        C = dict(zip(df.index, colors))
        node_colors = [C[node] if node in C else 'Black' for node in T.nodes()]
        #checks that colormap matches nodes so legend doesn't have crap in it
        cmap = check_keys(cmap, df[colorcol].unique())
        #add legend for node colors
        p = plotting.make_legend(ax.figure, cmap, loc=legend_loc, title=colorcol,fontsize=9)

    else:
        node_colors = 'black'
    nx.draw_networkx(T, pos, node_color=node_colors, node_size=node_size,
                     font_size=font_size, with_labels=with_labels, ax=ax)
    if edge_labels == True:
        nx.draw_networkx_edge_labels(T, pos, edge_labels=labels, font_size=font_size*.7, ax=ax)

    if labelcol not in [None,'']:
        node_labels = {node:df.loc[node][labelcol] if node in df.index else '' for node in T.nodes()}
        #print (node_labels)
        nx.draw_networkx_labels(T, pos, labels=node_labels, font_size=font_size,
                 horizontalalignment='left',verticalalignment='bottom')
    ax.axis('off')
    return T, pos

def check_keys(d, vals):
    """Remove dict keys not in vals"""

    new = d.copy()
    keys = list(d.keys())
    for key in keys:
        if not key in vals:
            new.pop(key, None)
    return new

def get_unique_snps(names, df, present=True):
    """Get snps unique to one or more samples from a SNP matrix.
    Args:
        name: name of sample(s)
        df: snp matrix from app.get_aa_snp_matrix(csq)
        present: whether snp should be present/absent
    returns:
        dataframe
    """

    if type(names) is str:
        names=[names]
    insamp = df[names]
    other = df.loc[:, ~df.columns.isin(names)]
    if present == True:
        u = other[other.sum(1)==0]
        u = df.loc[u.index]
    else:
        u = other[other.sum(1)==len(other.columns)]
        #sns.clustermap(df.loc[u.index])
        u = insamp.loc[u.index]
        u = u[u.sum(1)==0]
    return u

def get_fasta_length(filename):
    """Get length of reference sequence"""

    from pyfaidx import Fasta
    refseq = Fasta(filename)
    key = list(refseq.keys())[0]
    l = len(refseq[key])
    return l

def get_chrom(filename):
    """Get first chromosome name from fasta file"""

    rec = list(SeqIO.parse(filename, 'fasta'))[0]
    return rec.id

def get_fastq_info(filename):

    rl = get_fastq_read_lengths(filename)
    #print (filename,rl.mean())
    if rl.mean() == np.nan:
        return
    return int(rl.mean())

def get_fastq_num_reads(filename):
    """Return fastq mean number of reads"""

    cmd = 'expr $(zcat "%s" | wc -l) / 4' %filename
    tmp = subprocess.check_output(cmd, shell=True)
    l = int(tmp)
    return l

def get_fastq_read_lengths(filename, size=2000):
    """Return fastq read lengths"""

    df = fastq_to_dataframe(filename, size=size)
    #name = os.path.basename(filename).split('.')[0]
    return df.length

def get_file_size(filename):
    """Get a file size in MB"""

    stats1 = os.stat(filename)
    return round(stats1.st_size / (1024 * 1024),2)

def mafft_alignment(filename, outfile):
    """
    Align 2 sequences with mafft.
    Args:
        filename: fasta file with sequences to align
        seqs: sequences as list if no filename
    Returns:
        alignment object
    """

    cmd = f'mafft --auto {filename} > {outfile}'
    subprocess.check_output(cmd, shell=True)
    aln = AlignIO.read(outfile, 'fasta')
    return aln

def clustal_alignment(filename=None, seqs=None, command="clustalw"):
    """
    Align 2 sequences with clustal.
    Args:
        filename: fasta file with sequences to align
        seqs: sequences as list if no filename
        command: name of command, default clustalw
    Returns:
        alignment object
    """

    if filename == None:
        filename = 'temp.faa'
        SeqIO.write(seqs, filename, "fasta")
    name = os.path.splitext(filename)[0]
    from Bio.Align.Applications import ClustalwCommandline
    cline = ClustalwCommandline(command, infile=filename)
    stdout, stderr = cline()
    align = AlignIO.read(name+'.aln', 'clustal')
    return align

def seqrecords_from_alignment(alignment):
    """
    Remove gaps from sequences in a MultipleSeqAlignment
     object and return SeqRecord objects.
    """

    original_seqrecords = []
    for record in alignment:
        # Remove gaps from the sequence
        ungapped_seq = str(record.seq).replace('-', '')
        # Create a new SeqRecord with the ungapped sequence
        ungapped_record = SeqRecord(Seq(ungapped_seq), id=record.id, name=record.name, description=record.description)
        original_seqrecords.append(ungapped_record)

    return original_seqrecords

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
    #res['perc_qcov'] = res.length/len(res.qseq)*100
    return res

def local_blast(database, query, output=None, maxseqs=100, evalue=0.001,
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

def remote_blast(db, query, maxseqs=50, evalue=0.001, **kwargs):
    """Remote blastp.
    Args:
        query: fasta file with sequence to blast
        db: database to use - nr, refseq_protein, pdb, swissprot
    """
    from Bio.Blast.Applications import NcbiblastpCommandline
    output = os.path.splitext(query)[0]+'_blast.txt'
    outfmt = '"6 qseqid sseqid qseq sseq pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore stitle"'
    cline = NcbiblastpCommandline(query=query, db=db, max_target_seqs=maxseqs, outfmt=outfmt,
                                  evalue=evalue, out=output, remote=True, **kwargs)
    stdout, stderr = cline()
    return

def blast_fasta(database, filename, pident=90, **kwargs):
    """
    Args:
        pident: percent identity cutoff
    """

    remotedbs = ['nr','refseq_protein','pdb','swissprot']
    if database in remotedbs:
        remote_blast(database, filename, **kwargs)
    else:
        outfile = local_blast(database, filename, **kwargs)
    df = get_blast_results(filename=outfile)
    #filter with cutoffs
    df = df[(df.pident>=pident)]# & (df.perc_cov>=perc_cov)]
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
    df = blast_fasta(database, 'tempseq.fa', **kwargs)
    return df

def run_prokka(filename, path='.', trusted='', threads=4):
    """Run prokka on a fasta file.
        Args:
            outdir: results folder
            trusted: trusted proteins (gbk or fasta protein file)
    """

    name = os.path.splitext(os.path.basename(filename))[0]
    outdir = os.path.join(path, name)
    if trusted != '' and os.path.exists(trusted):
        cmd = f'/local/prokka/bin/prokka --proteins {trusted} \
            --outdir {outdir} {filename} --prefix {name} --centre X --compliant --cpus {threads}'
    else:
        cmd = f'/local/prokka/bin/prokka \
            --outdir {outdir} {filename} --prefix {name} --centre X --compliant --cpus {threads}'
    if not os.path.exists(outdir):
        print (cmd)
        subprocess.check_output(cmd, shell=True)
    else:
        print ('folder exists')
    return

def run_kraken(file1, file2='', dbname='STANDARD16', threads=4, outfile='krakenreport.txt'):
    """
    Run kraken2 on single/paired end fastq files. Requires that the dbname is set in
    your environment.
    """

    os.environ['KRAKEN2_DB_PATH'] = '/local/kraken2'
    if isinstance(file2, str) and os.path.exists(file2):
        #paired end
        cmd = 'kraken2 -db {db} --report {o} --threads {t} --confidence 0.1 --paired {f1} {f2} > kraken.out'\
            .format(f1=file1,f2=file2,t=threads,db=dbname,o=outfile)
    else:
        cmd = 'kraken2 -db {db} --report {o} --threads {t} --confidence 0.1 {f1} > kraken.out'\
            .format(f1=file1,t=threads,db=dbname,o=outfile)
    print (cmd)
    subprocess.check_output(cmd, shell=True)
    rep=pd.read_csv('krakenreport.txt',sep='\t',names=['perc_frag','n_frags_root','n_frags','rank_code','taxid','name'])
    rep['name'] = rep.name.str.lstrip()
    return rep

def kraken_run_samples(samples, dbname='STANDARD16', threads=4, found=pd.DataFrame({'sample':[]})):
    """
    Run a set of samples with kraken.
    """

    res=[]
    if not 'filename2' in samples.columns:
        samples['filename2'] = pd.Series()
    for i,r in samples.iterrows():
        name = r['sample']
        print (name)
        if name in list(found['sample']):
            continue
        try:
            rep = run_kraken(r.filename1,r.filename2, threads=threads, dbname=dbname)
        except Exception as e:
            print (e)
            pass
        rep = rep[rep.perc_frag>=0.1]
        rep['sample'] = name
        res.append(rep)

    if len(res)>0:
        df=pd.concat(res)
        found = pd.concat([found,df])
    return found

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
    df[key] = df[key].str.replace('|','_',regex=False)
    df['length'] = df.sequence.str.len()
    return df

def fastq_random_seqs(filename, size=50):
    """Random sequences from fastq file. Requires pyfastx.
    Creates a fastq index which will be a large file.
    """
    #see https://pythonforbiologists.com/randomly-sampling-reads-from-a-fastq-file.html
    import pyfastx
    fq = pyfastx.Fastq(filename, build_index=True)
    pos = np.random.randint(1,len(fq),size)
    seqs = [SeqRecord(Seq(fq[i].seq),id=str(fq[i].name)) for i in pos]
    return seqs

def fastq_to_rec(filename, size=50):
    """Get reads from a fastq file
        Args:
            size: limit
        Returns: biopython seqrecords
    """

    ext = os.path.splitext(filename)[1]
    if ext=='.gz':
        fastq_parser = SeqIO.parse(gzopen(filename, "rt"), "fastq")
    else:
        fastq_parser = SeqIO.parse(open(filename, "r"), "fastq")
    res=[]
    i=0
    for fastq_rec in fastq_parser:
        i+=1
        if i>size:
            break
        res.append(fastq_rec)
    return res

def fastq_to_dataframe(filename, size=5000):
    """Convert fastq to dataframe.
        size: limit to the first reads of total size, use None to get all reads
        Returns: dataframe with reads
    """

    if size==None:
        size=1e7
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
        keys = []
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
            keys += [i for i in quals.keys() if i not in keys]
        quals = keys+['id','start','end','strand','feat_type','sequence']
        df = pd.DataFrame(allfeat,columns=quals)
        if 'sequence' in df.keys():
            df['length'] = df.sequence.astype('str').str.len()
        res.append(df)
    res = pd.concat(res)
    if cds == True:
        res = res[res.feat_type=='CDS']
    if len(df) == 0:
        print ('ERROR: returned empty data, check that the file contains protein sequences '\
               'in the translation qualifier of each protein feature.' )
    return res

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

def trim_reads(filename1, filename2, outpath, quality=20,
                method='cutadapt', threads=4):
    """Trim adapters using cutadapt"""

    outfile1 = os.path.join(outpath, os.path.basename(filename1))
    outfile2 = os.path.join(outpath, os.path.basename(filename2))
    if method == 'cutadapt':
        cmd = 'cutadapt -m 80 -O 5 -q {q} -j {t} -o {o} -p {p}  {f1} {f2} '\
                .format(f1=filename1,f2=filename2,o=outfile1,p=outfile2,t=threads,q=quality)
        print (cmd)
        result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    return outfile1, outfile2

'''def get_subsample_reads(filename, outpath, reads=10000):
    """
    Sub-sample a fastq file with first n reads.
    Args:
        filename: input fastq.gz file
        outpath: output directory to save new file
        reads: how many reads to sample from start
    """

    lines = reads*4
    #print (lines % 4)
    name = os.path.basename(filename)
    #print (name)
    out = os.path.join(outpath, name)
    cmd = 'zcat {f} | head -n {r} | gzip > {o}'.format(f=filename,r=lines,o=out)
    print (cmd)
    subprocess.check_output(cmd, shell=True)
    return'''

def subset_reads(filename, path='/tmp', numreads=10000, overwrite=False):
    """Subset of reads"""

    from Bio import SeqIO, bgzf
    from gzip import open as gzopen
    new = os.path.join(path,os.path.basename(filename))
    if os.path.exists(new) and overwrite == False:
        return new
    print ('subsetting %s reads' %numreads)
    recs = SeqIO.parse(gzopen(filename, "rt"), "fastq")
    i=0
    with bgzf.BgzfWriter(new, "wb") as outgz:
        for rec in recs:
            if i>numreads:
                break
            SeqIO.write(sequences=rec, handle=outgz, format="fastq")
            i+=1
    return new

def get_vcf_samples(filename):
    """Get list of samples in a vcf/bcf"""

    from io import StringIO
    bcftoolscmd = get_cmd('bcftools')
    cmd = '{bc} query -l {f}'.format(bc=bcftoolscmd,f=filename)
    tmp = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    c = pd.read_csv(StringIO(tmp),sep='\t',names=['name'])
    return c.name

def vcf_to_dataframe(vcf_file, limit=1e5):
    """
    Convert a multi sample bcf/vcf to dataframe. Records each samples FORMAT fields.
    May use a lot of memory for large vcfs!
    Args:
        vcf_file: input multi sample vcf
        limit: limit to first n lines, avoids loading huge files
    Returns: pandas DataFrame
    """

    import vcf
    ext = os.path.splitext(vcf_file)[1]
    if ext == '.gz':
        file = gzopen(vcf_file, "rt")
    elif ext == '.bcf':
        # Convert BCF to VCF first
        bcftoolscmd = get_cmd('bcftools')
        cmd = '{bc} view {f} -O v -o temp.vcf'.format(bc=bcftoolscmd,f=vcf_file)
        tmp = subprocess.check_output(cmd, shell=True, universal_newlines=True)
        vcf_file = 'temp.vcf'
        file = open(vcf_file)
    else:
        file = open(vcf_file)
    vcf_reader = vcf.Reader(file,'r')
    res=[]
    cols = ['sample','REF','ALT','mut','DP','ADF','ADR','AD','chrom','var_type',
            'sub_type','pos','start','end','QUAL']
    i=0
    for rec in vcf_reader:
        if i>limit:
            break
        x = [rec.CHROM, rec.var_type, rec.var_subtype, rec.POS, rec.start, rec.end, rec.QUAL]
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
        i+=1

    res = pd.DataFrame(res,columns=cols)
    res = res[~res.start.isnull()]
    res['pos'] = res.pos.astype(int)
    #res['end'] = res.end.astype(int)
    return res

def get_vcf_positions(vcf_file):
    """Extract all positions from vcf per sample into a dataframe.
      Used in site proxmity filtering."""

    import vcf
    from gzip import open as gzopen
    ext = os.path.splitext(vcf_file)[1]
    if ext == '.gz':
        file = gzopen(vcf_file, "rt")
    else:
        file = open(vcf_file)
    vcf_reader = vcf.Reader(file,'r')

    res=[]
    cols = ['sample','pos']
    for rec in vcf_reader:
        x = []
        for sample in rec.samples:
            row = [sample.sample, rec.POS]
            res.append(row)

    res = pd.DataFrame(res,columns=cols)
    return res

def bcftools_query(bcf_file, positions=[], field='AD', show_cmd=False):
    """Query a vcf/bcf file for specific field over given positions.
        Args:
            positions: list of sites you want to query
            field: FORMAT field to extract e.g. DP, AD
    """

    bcftoolscmd = get_cmd('bcftools')
    pstr = ''
    for p in positions:
        pstr+='POS=%s|' %p
    pstr =pstr[:-1]
    cmd = "{bc} query -f '%CHROM %POS %REF %ALT [%{f} ]\\n' -i '{p}' {b}".format(b=bcf_file,
            bc=bcftoolscmd,f=field,p=pstr)
    tmp = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    if show_cmd == True:
        print (cmd)
    from io import StringIO
    df = pd.read_csv(StringIO(tmp),sep=' ',header=None)
    df = df.rename(columns={0:'chrom',1:'pos',2:'ref',3:'alt'})
    return df

def bcftools_count_sites(vcf_file):
    """Count no. of sites in vcf file"""

    bcftoolscmd = get_cmd('bcftools')
    cmd = "{bc} query -f 'x' {b}".format(b=vcf_file,bc=bcftoolscmd)
    tmp = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    return len(tmp)

def bcftools_call(bcf_file, vcfout, show_cmd=False):
    """Run bcftools call"""

    bcftoolscmd = get_cmd('bcftools')
    cmd = '{bc} call --ploidy 1 -m -v -o {o} {raw}'.format(bc=bcftoolscmd,o=vcfout,raw=bcf_file)
    if show_cmd == True:
        print (cmd)
    tmp = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    return

def bcftools_filter(vcf_file, out=None, filters=''):
    """Count no. of sites in vcf file"""

    bcftoolscmd = get_cmd('bcftools')
    if out == None:
        cmd = '{bc} filter -i "{f}" {i}'.format(bc=bcftoolscmd,i=vcf_file,f=filters)
    else:
        cmd = '{bc} filter -i "{f}" -o {o} -O z {i}'.format(bc=bcftoolscmd,i=vcf_file,o=out,f=filters)
    print (cmd)
    tmp = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    return

def bcftools_merge(bcf_files, out, threads=4):
    """Merge bcf files"""

    bcftoolscmd = get_cmd('bcftools')
    cmd = '{bc} merge --threads {t} -m all -o {o} {b}'.format(b=bcf_files, o=out, bc=bcftoolscmd, t=threads)
    print (cmd)
    subprocess.check_output(cmd, shell=True)
    return

def bcftools_consensus(vcf_file, ref, out):
    """Get consensus sequences from output of mpileup
    see https://samtools.github.io/bcftools/howtos/consensus-sequence.html
     """

    bcftoolscmd = get_cmd('bcftools')
    cmd = '{bc} index {f}'.format(bc=bcftoolscmd,f=vcf_file)
    tmp = subprocess.check_output(cmd,shell=True)
    cmd = 'cat {r} | {bc} consensus {f} > {o}'.format(bc=bcftoolscmd,r=ref,o=out,f=vcf_file)
    print (cmd)
    tmp = subprocess.check_output(cmd,shell=True)
    print ('consensus sequences saved to %s' %out)
    return

def spades(file1, file2, path, outfile=None, threads=4, trusted=None, cmd='spades'):
    """Run spades on paired end reads"""

    if trusted != None:
        cmd = '%s --careful -t %s --trusted-contigs %s --pe1-1 %s --pe1-2 %s --careful -o %s' %(cmd,threads,trusted,file1,file2,path)
    else:
        cmd = '%s --careful -t %s --pe1-1 %s --pe1-2 %s --careful -o %s' %(cmd,threads,file1,file2,path)
    if not os.path.exists(path):
        print (cmd)
        subprocess.check_output(cmd, shell=True)
    if outfile != None:
        shutil.copy(os.path.join(path,'scaffolds.fasta'),outfile)
    return outfile

def get_snp_matrix(df):
    """SNP matrix from multi sample vcf dataframe"""

    df = df.drop_duplicates(['mut','sample'])
    x = df.set_index(['mut','sample'])['start'].unstack('sample')
    x[x.notna()] = 1
    x = x.fillna(0)
    return x

def plot_fastq_qualities(filename, ax=None, limit=20000,
                         title='per base sequence quality'):
    """Plot fastq qualities for illumina reads."""

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
    boxprops = dict(linestyle='-', linewidth=.5, color='black')
    df.plot(kind='box', ax=ax, grid=False, showfliers=False, boxprops=boxprops,
            color=dict(boxes='black',whiskers='black'))
    n = int(l/10)
    ax.set_xticks(np.arange(0, l, n))
    ax.set_xticklabels(np.arange(0, l, n))

    ax.set_xlabel('position(bp)')
    ax.set_xlim((0,l))
    ax.set_ylim((0,40))
    ax.set_title(title, fontsize=10)
    return

def normpdf(x, mean, sd):
    """Normal distribution function"""

    import math
    var = float(sd)**2
    denom = (2*math.pi*var)**.5
    num = math.exp(-(float(x)-float(mean))**2/(2*var))
    return num/denom

def get_fasta_gc(fasta_file, limit=1e4):
    """Find gc content across sequence file"""

    from Bio.SeqUtils import gc_fraction
    #seq = SeqIO.read(fasta_file, format='fasta')
    recs = list(SeqIO.parse(fasta_file, "fasta"))
    seq = recs[0]
    gc = gc_fraction(seq)*100
    return gc

def get_fastq_gc(fastq_file, limit=1e4):
    """Find gc content across sequence file"""

    from Bio.SeqUtils import gc_fraction
    df = fastq_to_dataframe(fastq_file, size=limit)
    gc = df.seq.apply(lambda x: gc_fraction(x)*100)
    return gc

def plot_fastq_gc_content(filename, ax=None, limit=50000):
    """Plot fastq gc conent"""

    import pylab as plt
    from Bio.SeqUtils import gc_fraction
    if ax == None:
        f,ax=plt.subplots(figsize=(12,5))
    plt.style.use('bmh')
    df = fastq_to_dataframe(filename, size=limit)
    gc = df.seq.apply(lambda x: gc_fraction(x)*100)
    gc.hist(ax=ax,bins=80,color='black',grid=False,histtype='step',lw=1)
    ax.set_xlim((0,100))
    x = np.arange(1,100,.1)
    meangc = round(gc.mean(),2)
    f = [normpdf(i, meangc, gc.std()) for i in x]
    ax2 = ax.twinx()
    ax2.plot(x,f)
    ax2.set_ylim(0,max(f))
    ax.set_title('GC content. Mean=%s' %meangc,size=12)
    return

def get_instrument(fastq_file):

    import gzip, re
    instrumentIDs = {"HWI-M[0-9]{4}" : "MiSeq",
            "HWUSI" : "Genome Analyzer IIx",
            "M[0-9]{5}" : "MiSeq",
            "HWI-C[0-9]{5}" : "HiSeq 1500",
            "C[0-9]{5}" : "HiSeq 1500",
            "HWI-D[0-9]{5}" : "HiSeq 2500",
            "D[0-9]{5}" : "HiSeq 2500",
            "J[0-9]{5}" : "HiSeq 3000",
            "K[0-9]{5}" : "HiSeq 3000/HiSeq 4000",
            "E[0-9]{5}" : "HiSeq X",
            "NB[0-9]{6}": "NextSeq",
            "NS[0-9]{6}" : "NextSeq",
            "MN[0-9]{5}" : "MiniSeq",
            "VH[0-9]{5}" : "NextSeq",
            "H[A-Z,0-9]{4}": "NovaSeq S4"}

    with gzip.open(fastq_file) as f:
        head = ""
        line = str(f.readline())
    #print (line)
    for pattern,iid in instrumentIDs.items():
        s = re.search(pattern, line)
        if s != None:
            #print (iid)
            return iid
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

def concatenate_fasta(input_files, output_file, max_length=1e7):
    """Join fasta files into one"""

    with open(output_file, 'w') as output_handle:
        for input_file in input_files:
            for record in SeqIO.parse(input_file, 'fasta'):
                if len(record.seq) > max_length:
                    print (record.id, len(record.seq), 'too large')
                    continue
                SeqIO.write(record, output_handle, 'fasta')
    return

def core_alignment_from_vcf(vcf_file, uninformative_sites=False, missing=False, omit=None, callback=None):
    """
    Get core SNP site calls as sequences from a multi sample vcf file.
    Args:
        vcf_file: multi-sample vcf (e.g. produced by app.variant_calling)
        uninformative_sites: whether to include uninformative sites
        missing: whether to include sites with one or more missing samples (ie. no coverage)
        omit: list of samples to exclude if required
    """

    import vcf
    from collections import defaultdict
    vcf_reader = vcf.Reader(open(vcf_file, 'rb'))
    print ('uninformative_sites', uninformative_sites)
    def default():
        return []
    result = defaultdict(default)
    sites = []
    result['ref'] = []
    missing_sites = []
    uninf_sites = []
    for record in vcf_reader:
        S = {sample.sample: sample.gt_bases for sample in record.samples}
        #print (record.POS,S)
        if omit != None:
            for o in omit:
                del S[o]
        #if any missing samples at the site we don't add
        if None in S.values():
            #print (S)
            missing_sites.append(record.POS)
            if missing == False:
                continue
        #ignore uninformative sites
        if uninformative_sites == False and len(S)>1:
            u = set(S.values())
            if len(u) == 1:
                uninf_sites.append(record.POS)
                continue
        result['ref'].append(record.REF)
        #get bases over all samples
        for name in S:
            val = S[name]
            if val == None:
                val = 'N'
            result[name].append(val)
        sites.append(record.POS)

    sites = sorted(list(set(sites)))
    print ('found %s sites for core snps' %len(sites))
    print ('%s sites with at least one missing sample' %len(missing_sites))
    if uninformative_sites == False:
        print ('%s uninformative sites' %len(uninf_sites))
    if len(sites)==0:
        print ('no sites found may mean:\n'
         '- samples are identical\n'
         '- one sample is too different due to poor coverage/depth\n'
         'this may be caused by contamination or low quality data')
    recs = []
    for sample in result:
        seq = ''.join(result[sample])
        seqrec = SeqRecord(Seq(seq),id=sample)
        recs.append(seqrec)
        #print (len(seqrec))

    smat = pd.DataFrame(result)
    smat.index = sites
    smat.index.rename('pos', inplace=True)
    smat = smat.sort_index()
    return recs, smat

def samtools_flagstat(filename):
    """Parse samtools flagstat output into dictionary"""

    samtoolscmd = get_cmd('samtools')
    cmd = '{s} flagstat {f}'.format(f=filename,s=samtoolscmd)
    tmp = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    x = tmp.split('\n')
    x = [int(i.split('+')[0]) for i in x[:-1]]
    #print (x)
    cols = ['total','primary','secondary','supplementary','duplicates','primary duplicates',
            'mapped','primary mapped',
            'paired','read1','read2','properly paired','with itself','singletons']
    d = {}
    for c,v in zip(cols,x):
        #print (c,v)
        d[c] = v

    #d['perc_mapped'] = round(d['mapped']/d['total']*100,2)
    return d

def samtools_tview(bam_file, chrom, pos, width=200, ref='', display='T'):
    """View bam alignment with samtools"""

    samtoolscmd = get_cmd('samtools')
    os.environ["COLUMNS"] = str(width)
    cmd = '{sc} tview {b} -p {c}:{p} -d {d} {r}'\
            .format(b=bam_file,c=chrom,p=pos,d=display,r=ref,sc=samtoolscmd)
    print (cmd)
    tmp = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    return tmp

def samtools_coverage(bam_file):
    """Get coverage/depth stats from bam file"""

    samtoolscmd = get_cmd('samtools')
    cmd = '{sc} coverage {b}'.format(b=bam_file,sc=samtoolscmd)
    tmp = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    from io import StringIO
    c = pd.read_csv(StringIO(tmp),sep='\t',header=0).round(2)
    c = c.iloc[0]
    return c

def samtools_depth(bam_file, chrom=None, start=None, end=None):
    """Get depth from bam file"""

    samtoolscmd = get_cmd('samtools')
    if start != None and end != None:
        if chrom == None:
            #try to get chromosome name
            chrom = get_chrom_names(bam_file)
        cmd = '{sc} depth -r {c}:{s}-{e} {b}'.format(b=bam_file,c=chrom,s=start,e=end,sc=samtoolscmd)
    else:
        cmd = '{sc} depth {b}'.format(b=bam_file,sc=samtoolscmd)

    tmp = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    c = pd.read_csv(StringIO(tmp),sep='\t',names=['chr','pos','depth'])
    return c

def get_mean_depth(bam_file, chrom=None, start=None, end=None, how='mean'):
    """Get mean depth from bam file"""

    c = samtools_depth(bam_file, chrom, start, end)
    if how == 'mean':
        return c.depth.mean().round(2)
    else:
        return c.depth.sum().round(2)

def get_chrom_names(bam_file):
    """Get chromosome names from bam file"""

    samtoolscmd = get_cmd('samtools')
    cmd = '{sc} idxstats {b}'.format(b=bam_file,sc=samtoolscmd)
    tmp = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    df = pd.read_csv(StringIO(tmp),sep='\t',names=['chrom','start','end','xx'])
    return df.iloc[0].chrom

def fetch_sra_reads(df, path, key='Run', test=False, sample_size=None, omit=[]):
    """Download a set of reads from SRA using dataframe with runs"""

    if sample_size != None:
        df = df.sample(sample_size)
    for i,r in df.iterrows():
        name = r[key]
        if name in omit: continue
        files = glob.glob(os.path.join(path,name+'*'))
        #if files exist we don't download
        if len(files) == 0:
            cmd = 'fastq-dump --split-files {n} -O {o}'.format(n=r[key],o=path)
            print (cmd)
            if test == False:
                subprocess.check_output(cmd,shell=True)
    return

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

    def clean_location(location):
        """Remove < or > from feature location coordinates."""
        start = location.start
        end = location.end
        # Ensure start and end are integers and remove special characters
        start = int(str(start).replace("<", "").replace(">", ""))
        end = int(str(end).replace("<", "").replace(">", ""))
        return FeatureLocation(start, end, strand=location.strand)

    l=1
    prefix = create_locus_tag(in_file)
    recs = SeqIO.parse(in_file,format='gb')
    for record in recs:
        #make a copy of the record as we will be changing it during the loop
        new = copy(record)
        new.features = []
        #loop over all features
        for i in range(0,len(record.features)):
            feat = record.features[i]
            feat.location = clean_location(feat.location)
            q = feat.qualifiers
            #remove some unecessary qualifiers
            for label in ['note','translation','product','experiment']:
                if label in q:
                    del q[label]
            if(feat.type == "CDS"):
                #use the CDS feature to create the new lines
                if 'locus_tag' in q:
                    tag = q['locus_tag'][0]
                elif 'gene' in q:
                    tag = prefix+'_'+q['gene'][0]
                    q['locus_tag'] = tag
                else:
                    tag = q['locus_tag'] = i
                q['ID'] = 'CDS:%s' %tag
                q['Parent'] = 'transcript:%s' %tag
                q['biotype'] = 'protein_coding'
                #create mRNA feature
                m = SeqFeature(feat.location,type='mRNA')#,strand=feat.strand)
                q = m.qualifiers
                q['ID'] = 'transcript:%s' %tag
                q['Parent'] = 'gene:%s' %tag
                q['biotype'] = 'protein_coding'
                new.features.append(m)
            elif(record.features[i].type == "gene"):
                #edit the gene feature
                q=feat.qualifiers
                #print(feat.location.start)
                #print (feat)
                if 'locus_tag' in q:
                    tag = q['locus_tag'][0]
                else:
                    tag = prefix+'_'+q['gene'][0]
                    q['locus_tag'] = tag
                q['ID'] = 'gene:%s' %tag
                q['biotype'] = 'protein_coding'
                if 'gene' in q:
                    q['Name'] = q['gene']
            new.features.append(feat)
        #write the new features to a GFF
        GFF.write([new], out_handle)
        return

def create_locus_tag(filename):
    """Create a genbank style locus tag"""

    name = os.path.basename(filename)
    x = hashlib.md5(name.encode('utf-8')).hexdigest()
    x = ''.join(i for i in x if not i.isdigit())
    x = x.upper()[:7]
    return x

def get_spoligotype(filename, reads_limit=3000000, threshold=2, threads=4):
    """
    Get spoligotype from WGS reads. Blasts to known DR spacers.
    Args:
        reads_limit: limit of subset of reads to blast
        threshold: minimum number of hits to a spacer to consider it present
    """

    ref = os.path.join(datadir, 'dr_spacers.fa')
    #convert reads to fasta
    try:
        fastq_to_fasta(filename, 'temp.fa', reads_limit)
    except Exception as e:
        print(e)
        return
    #make blast db from reads
    make_blast_database('temp.fa')
    #blast spacers to db
    bl = blast_fasta('temp.fa', ref, evalue=0.1, threads=threads,
                           maxseqs=reads_limit, show_cmd=False)
    bl=bl[(bl.qcovs>90) & (bl.mismatch<=threshold)]
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

def get_spoligotypes(samples, spo=None):
    """Get spoligotypes for multiple M.bovis strains"""

    if spo is not None:
        done=list(spo['sample'])
    else:
        done=[]
    samples = samples.drop_duplicates('sample')
    res=[]
    for i,r in samples.iterrows():
        f=r.filename1
        samp=r['sample']
        if samp in done:
            continue
        b = get_spoligotype(f)
        sb = get_sb_number(b)
        print (r['sample'], sb, b)
        res.append([r['sample'],sb,b])

    res = pd.DataFrame(res,columns=['sample','SB','code'])
    return res

def make_mask_file(gb_file, outfile):
    """
    Make mask bed file from genbank annoation by finding repeat regions.
    """

    g = genbank_to_dataframe(gb_file)
    g['gene'] = g.gene.fillna('')
    g = g[((g.gene.str.contains('PE') | g.gene.str.contains('PPE')) & (g.feat_type=='CDS')) | (g.feat_type=='repeat_region')]
    if len(g)==0:
        print ('no regions found to mask')
        return
    lines = []
    #get chromosome/genome name (assumes only one)
    chrom = g.iloc[0].id
    for i,r in g.iterrows():
        l = '{c} \t {s}\t{e}\t{g}\t{t}'.format(c=chrom,s=r.start,e=r.end,g=r.gene,t=r.locus_tag)
        #print (l)
        lines.append(l)
    print ('found %s regions' %len(lines))
    with open(outfile, 'w') as fp:
        for l in lines:
            fp.write("%s\n" %l)
    return lines

def compare_results(c1, c2, sample=None):
    """
    Compare two runs of snipgenie.
    """

    x = c1[~c1.pos.isin(c2.pos)]
    y = c2[~c2.pos.isin(c1.pos)]
    print ('%s/%s sites not in second:' %(len(x),len(c1)))
    print (x)
    print ('-------------------------')
    print ('%s/%s sites not in first:' %(len(y),len(c2)))
    print (y)
    return x,y

def calculate_N50(list_of_lengths):
    """Calculate N50 for a sequence of numbers. For assemblies.
    Args:
        list_of_lengths (list): List of numbers.
    Returns:
        float: N50 value.
    """
    tmp = []
    for tmp_number in set(list_of_lengths):
            tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
    tmp.sort()
    if (len(tmp) % 2) == 0:
        median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
    else:
        median = tmp[int(len(tmp) / 2)]
    return median

def get_sequence_lengths(fasta_file):
    """Get sequence lengths from fasta"""

    lengths = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        lengths.append(len(record.seq))
    return lengths
