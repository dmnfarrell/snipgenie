"""
    MLST methods for M.bovis.
    Created May 2021
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
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from . import tools, app

home = os.path.expanduser("~")
config_path = os.path.join(home,'.config','snipgenie')
module_path = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(module_path, 'data')
mlstdir = os.path.join(module_path, 'mlst')

mbovis_scheme = pd.read_csv(os.path.join(mlstdir,'mbovis_scheme.csv'))
mbovis_db = os.path.join(mlstdir,'mbovis_db.csv.gz')
sample_profiles = os.path.join(mlstdir,'sample_profiles.csv')
ref_proteins = os.path.join(datadir,'Mbovis_AF212297_proteins.fa')

schemes = {'Mbovis-AF212297': mbovis_scheme}

def get_samples_vcf(vcf_file):
    """Get all samples in a vcf file. Returns a list."""

    cmd = 'bcftools query -l %s' %vcf_file
    tmp = subprocess.check_output(cmd, shell=True)
    names = tmp.decode().split('\n')[:-1]
    #print (names)
    return names

def get_nucleotide_sequences(gb_file,out_file,idkey='locus_tag'):
    """Extract nucleotide sequences for all features in genbank file"""

    recs = SeqIO.to_dict(SeqIO.parse(gb_file,'genbank'))
    chroms = list(recs.keys())
    result = []
    for chrom in chroms:
        rec = recs[chrom]
        for f in rec.features[1:]:
            q=f.qualifiers
            if f.type != 'CDS':
                continue
            seq = rec.seq[f.location.start:f.location.end]
            try:
                new = SeqRecord(seq,id=q[idkey][0])
                result.append(new)
            except:
                #print (q)
                pass
    SeqIO.write(result,out_file,format='fasta')
    return result

def bcftools_consensus(vcf_file, sample, out_file='consensus.fa'):
    """Get consensus sequence from a vcf file"""

    cmd='bcftools index -f %s' %vcf_file
    subprocess.check_output(cmd, shell=True)
    cmd='bcftools consensus -f {r} -s {s} {v} > {o}'.format(r=app.mbovis_genome,
                                                v=vcf_file,s=sample,o=out_file)
    print (cmd)
    subprocess.check_output(cmd, shell=True)
    return

def make_reference_proteins():
    """Make reference protein sequences from m.bovis
       used for initial setup of proteins file which is used as
      'trusted source' in annotation process
     """

    prots = pg.tools.genbank_to_dataframe(app.mbovis_gb,cds=True)
    prots = prots.fillna('')
    prots = prots.dropna(subset=['locus_tag'])
    ref_proteins = 'Mbovis_AF212297_proteins.fa'
    #get prokka type header for using in annotation
    prots['header'] = prots.apply(lambda x: '~~~'.join([x.locus_tag,x.gene,x['product'],'none']),1)
    pg.tools.dataframe_to_fasta(prots,idkey='header',outfile=ref_proteins)
    print (len(prots))
    return prots

def find_alleles(fastafile):
    """Find alleles by simple matches to the reference table of known sequences.
    Checks if an allele already exists and if not assigns a new number.
    Returns:
        dataframe with allele number for each gene
        dataframe with new alleles to add to db
    """

    db = pd.read_csv(mbovis_db)
    names = db.name.unique()
    df = tools.fasta_to_dataframe(fastafile).reset_index()
    print (len(df))
    result=[]
    new=[]
    missed=0
    for name in names:
        #print (name)
        s = db[db.name==name]
        gene = df[df.name==name]
        #print (gene)
        if len(gene)==0:
            print ('missing gene:',name)
            #missing gene in target
            result.append((name,0))
            missed+=1
            continue
        target = gene.iloc[0].sequence
        found = s[s.sequence==target]
        #print (target,found)
        if len(found)>0:
            found = found.iloc[0]
            result.append((name,found.allele))
        else:
            #assign new allele
            newallele = s.allele.max()+1
            result.append((name,newallele))
            new.append([name,newallele,target])
    prof = pd.DataFrame(result,columns=['name','allele'])
    prof['allele'] = prof.allele.astype(int)
    #new additions
    new = pd.DataFrame(new,columns=['name','allele','sequence'])
    print ('missed', missed)
    return prof, new

def update_mlst_db(new):
    """Update the database of MLST profiles"""

    db = pd.read_csv(mbovis_db)
    #check if any new alleles to be added first
    found = new[new.sequence.isin(db.sequence)]
    if len(found) > 0:
        db = pd.concat([db,new])
    db.to_csv(mbovis_db, index=False, compression='gzip')
    print ('added %s new alleles' %len(new))
    return

def type_sample(fastafile, annotfile, outfile, threads=4, overwrite=False, update=True):
    """
    Type a single sample using wgMLST method. Requires pathogenie for
    annotation steps.
    Args:
        fastafile: fasta file to type from assembly or other
        annotfile: annotation file
        outfile: output file for annotations
        update: whether to update DB, default True
    Returns:
        dataframe containing the MLST profile
    """

    import pathogenie
    if overwrite == True or not os.path.exists(annotfile):
        #annotate
        featdf,recs = pathogenie.run_annotation(fastafile, threads=threads,
                                        kingdom='bacteria', trusted=ref_proteins)
        #get nucl sequences from annotation
        SeqIO.write(recs,annotfile,'genbank')
    if overwrite == True or not os.path.exists(outfile):
        get_nucleotide_sequences(annotfile, outfile, idkey='protein_id')

    #find alleles
    profile,new = find_alleles(outfile)
    #update db
    if update == True:
        update_mlst_db(new)
    return profile

def diff_profiles(s1, s2):
    return sum(1 for a, b in zip(s1, s2) if a != b)

def get_profile_string(df):
    return ';'.join(df.allele.astype(str))

'''def dist_matrix(profiles):
    """Distance matrix of a set of profiles"""

    dist=[]
    for s in profiles:
        x = profiles[s]
        row=[]
        for s in profiles:
            d = diff_profiles(x,profiles[s])
            row.append(d)
        dist.append(row)
    D = pd.DataFrame(dist,columns=profiles.keys(),index=profiles.keys())
    return D'''

def dist_matrix(profiles):
    """
    Distance matrix from a set of allele numbers.
    profiles: dataframe of allele profiles
    returns:
        dataframe
    """

    dist=[]
    for s,r in profiles.iterrows():
        x = list(r)
        row=[]
        for i,r2 in profiles.iterrows():
            d = diff_profiles(x,list(r2))
            row.append(d)
        dist.append(row)
    
    D = pd.DataFrame(dist,columns=profiles.index,index=profiles.index)
    return D

def tree_from_distmatrix(D):
    """tree from distance matrix"""

    from skbio import DistanceMatrix
    from skbio.tree import nj
    ids = list(D.index)
    dm = DistanceMatrix(D.values, ids)
    tree = nj(dm)
    #print(tree.ascii_art())
    return tree

'''def run_samples(vcf_file, outdir, names=None, omit=[], **kwargs):
    """Run all the samples in a vcf file.
    Args:
        vcf_file: multi sample variant file from previous calling
        outdir: folder for writing intermediate files
    Returns:
        dict of mst profiles
    """

    profs = {}
    if names == None:
        names = get_samples_vcf(vcf_file)

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    for sample in names:
        print (sample)
        if sample in omit:
            continue
        #get consensus sequences
        consfile = os.path.join(outdir, '%s_consensus.fa' %sample)
        if not os.path.exists(consfile):
            bcftools_consensus(vcf_file, sample, consfile)
        #seq = SeqIO.read('consensus.fa','fasta')
        #print (len(seq))
        annotfile = os.path.join(outdir, '%s.gb' %sample)
        outfile = os.path.join(outdir, '%s.fa' %sample)
        profile = type_sample(consfile, annotfile, outfile, **kwargs)
        profs[sample] = list(profile.allele)
    #convert to dataframe
    profs = pd.DataFrame(profs).T
    return profs'''

def run_samples(samples, outdir, names=None, omit=[], threads=4, **kwargs):
    """Run all the samples in a vcf file.
    Args:
        samples (dataframe): table of samples with fastq file names
        outdir: folder for writing intermediate files
    Returns:
        dict of mst profiles
    """

    profs = {}

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    
    if names!=None:
        samples = samples[samples['sample'].isin(names)]
    for i,r in samples.iterrows():    
        name=r['sample']        
        if name in omit:
            continue
        print (name)
        #run assembly here       
        assembly = os.path.join(outdir, '%s.fa' %name)
        print (assembly)
        tools.spades(r.filename1, r.filename2, os.path.join(outdir,name), outfile=assembly, threads=12)

        annotfile = os.path.join(outdir, '%s.gb' %name)
        outfile = os.path.join(outdir, '%s_proteins.fa' %name)
        profile = type_sample(assembly, annotfile, outfile, **kwargs)
        profs[name] = list(profile.allele)
    #convert to dataframe
    profs = pd.DataFrame(profs).T
    return profs

def main():
    "Run the application"

    import sys, os
    from argparse import ArgumentParser
    parser = ArgumentParser(description='SNiPgenie wgMLST tool. https://github.com/dmnfarrell/snipgenie')
    parser.add_argument("-i", "--input", dest="vcf_file", default=None,
                        help="input multi-sample vcf from variant calling", metavar="FILE")
    parser.add_argument("-o", "--outdir", dest="outdir",
                        help="Folder to output intermediate results", metavar="FILE")
    #parser.add_argument("-S", "--species", dest="species", default=None,
    #                    help="set the species")
    parser.add_argument("-t", "--threads", dest="threads", default=None,
                        help="cpu threads to use")
    parser.add_argument("-x", "--test", dest="test",  action="store_true",
                        default=False, help="Test run")

    args = vars(parser.parse_args())
    if args['test'] == True:
        test()
    elif args['vcf_file'] != None:
        profs = run_samples(args['vcf_file'], outdir=args['outdir'])
        D = dist_matrix(profs)
        D.to_csv('dist_mlst.csv',index=False)
        tree = tree_from_distmatrix(D)
        print()
        print(tree.ascii_art())
        tree.write('mlst.newick', 'newick')
    else:
        print('provide an input file using -i, use -h for help')

if __name__ == '__main__':
    main()
