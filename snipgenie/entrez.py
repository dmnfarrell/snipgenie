"""
    Various methods for bacterial genomics.
    Created Nov 2023
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
import urllib
from io import StringIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Phylo, AlignIO, Align
from Bio import Entrez
import numpy as np
import pandas as pd
import pylab as plt
import matplotlib as mpl

def get_assembly_summary(id):
    """Get esummary for an entrez id"""

    from Bio import Entrez
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record

def find_assemblies(term):
    """Find all assemblies for a search term, returns ids (ui list)
    Args:
        term: search term, usually organism name
    """
    #provide your own mail here
    Entrez.email = "A.N.Other@example.com"
    handle = Entrez.esearch(db="assembly", term=term, retmax='5000')
    record = Entrez.read(handle)
    ids = record['IdList']
    print (f'found {len(ids)} ids')
    return ids

def get_assembly_report(ids):
    """Download genbank assembly meta data for a set of ids, including
    links to download.
    Args:
        ids: ui list
    """

    from tqdm import tqdm
    Entrez.email = "A.N.Other@example.com"
    links = []
    result = []
    for id in tqdm(ids):
        row = {'id':id}
        #get summary
        rec = get_assembly_summary(id)
        #print (id)
        asm_summ = rec['DocumentSummarySet']['DocumentSummary'][0]
        fields = ['AssemblyAccession','BioSampleAccn','BioSampleId','SubmitterOrganization']
        for key in fields:
            row[key] = asm_summ[key]
        row['GenbankAccession'] = asm_summ['Synonym']['Genbank']

        #biosample info is a separate request using the BioSampleId
        handle = Entrez.esummary(db="biosample", id=asm_summ['BioSampleId'], report="full")
        rec2 = Entrez.read(handle)
        sampledata = rec2['DocumentSummarySet']['DocumentSummary'][0]['SampleData']
        #parse xml in sampledata
        from bs4 import BeautifulSoup
        soup = BeautifulSoup(sampledata)
        all_attr = soup.findAll('attribute')
        for attr in all_attr:
            #print (attr,attr['attribute_name'],attr.text)
            row[attr['attribute_name']] = attr.text
        #get url
        #url = asm_summ['FtpPath_RefSeq']
        url = asm_summ['FtpPath_GenBank']
        #print (row)
        if url != '':
            label = os.path.basename(url)
            #get the fasta link - change this to get other formats
            link = os.path.join(url,label+'.fna.gz')
            row['link'] = link
        result.append(row)     
    result = pd.DataFrame(result)  
    return result

def download_links(df, path):
    """download links from assemblue table"""

    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
    for i,r in df.iterrows():
        label = r.AssemblyAccession
        print (r.link)
        urllib.request.urlretrieve(r.link, os.path.join(path, f'{label}.fna.gz'))
