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
from . import tools

home = os.path.expanduser("~")
config_path = os.path.join(home,'.config','snipgenie')
module_path = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(module_path, 'data')

mlst = pd.read_csv(os.path.join(datadir,'mlst_scheme.csv'))
db = pd.read_csv(os.path.join(datadir,'mlst_db.csv.gz'))
ref_proteins = os.path.join(datadir,'Mbovis_AF212297_proteins.fa')

def get_samples_vcf(vcf_file):
    cmd = 'bcftools query -l %s' %vcf_file
    tmp = subprocess.check_output(cmd, shell=True)
    return tmp.decode().split('\n')

def get_consensus(vcf_file, sample, out_file):
    """Get consensus sequence from vcf"""

    cmd='bcftools index %s' %vcf_file
    #subprocess.check_output(cmd, shell=True)
    cmd='cat {r} | bcftools consensus -s {s} {v} > {o}'.format(r=snpg.app.mbovis_genome,v=vcf_file,s=sample,o=out_file)
    print (cmd)
    subprocess.check_output(cmd, shell=True)
    return

def get_profile_string(df):
    return ''.join(df.allele.astype(str))

def diff_profiles(s1, s2):
    return sum(1 for a, b in zip(list(s1), list(s2)) if a != b)

def find_alleles(fastafile):
    """Find allele by simple matches to the reference table of known sequences.
    Returns:
        dataframe with allele number for each gene
        dataframe with new alleles to add to db
    """

    db = pd.read_csv('mlst_db.csv')
    names = ref.name.unique()
    df = pg.tools.fasta_to_dataframe(fastafile)
    result=[]
    new=[]
    for name in names[:100]:
        #print (name)
        s = db[db.name==name]
        target = df[df.name==name].iloc[0].sequence
        found = s[s.sequence==target]
        #print (found)
        if len(found)>0:
            found = found.iloc[0]
            result.append((name,found.allele))
        else:
            #assign new allele
            newallele = s.allele.max()+1
            result.append((name,newallele))
            new.append([name,newallele,target])
    res = pd.DataFrame(result,columns=['name','allele'])
    new = pd.DataFrame(new,columns=['name','allele','sequence'])
    return res, new

def update_mlst_db(new):
    """Update the database of MLST profiles"""

    db = pd.read_csv('mlst_db.csv')
    db = pd.concat([db,new])
    db.to_csv('mlst_db.csv',index=False)
    print ('added %s new alleles' %len(new))
    return

def type_sample(sample, vcf_file, path, threads=4, overwrite=False):
    """Type a single sample using wgMLST.
    Args:
        sample: sample name, must be present in vcf file
        vcf_file: source vcf
        path: output folder for annotations
    Returns:
        dataframe of MLST profile
    """

    seqfile = 'consensus.fa'
    fastafile = os.path.join(path,'%s.fa' %sample)
    #use consensus for now
    get_consensus(vcf_file, sample, seqfile)
    if overwrite == True or not os.path.exists(fastafile):
        #annotate
        featdf,recs = pg.run_annotation(seqfile,
                                    threads=threads, kingdom='bacteria', trusted=ref_proteins)
        #get nucl sequences from annotation
        SeqIO.write(recs,'temp.gb','genbank')
        get_nucleotide_sequences('temp.gb',fastafile,idkey='protein_id')
    #find alleles
    res,new = find_alleles(fastafile)
    #print (res)
    #update db
    update_mlst_db(new)
    return res

def dist_matrix(profiles):
    """Distance matrix of a set of profiles"""

    dist=[]
    for s in profiles:
        x=profs[s]
        row=[]
        for s in profs:
            d = diff_profiles(x,profs[s])
            row.append(d)
        dist.append(row)
    D = pd.DataFrame(dist,columns=profs.keys(),index=profs.keys())
    return D
