"""
    Clustering methods for snipgenie.
    Created Feb 2023
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

import sys,os,subprocess,glob,shutil,re,random
import platform
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Phylo, AlignIO
import numpy as np
import pandas as pd
from  . import tools

snp200_cmap = {0:'darkgreen',1:'coral',2:'dodgerblue',3:'crimson',4:'lightgreen',
                5:'orange',6:'mediumpurple',7:'gold',8:'cyan',9:'beige',10:'pink',
                11:'brown',12:'gray',13:'blue',14:'green',15:'Aquamarine',
                16:'DarkCyan',17:'red',18:'IndianRed'}


'''def treecluster(filename, threshold=12):
    """Run treecluster methods on newick tree"""

    from treeswift import read_tree_newick
    infile = open(filename)
    for line in infile:
            if isinstance(line,bytes):
                l = line.decode().strip()
            else:
                l = line.strip()

    tree = read_tree_newick(l)
    clusters = min_clusters_threshold_max_clade(tree,threshold,0)
    cl = get_treecluster_labels(clusters)
    return cl.set_index('id')'''

def treecluster(filename, threshold=12, method='avg_clade'):
    """Run treecluster methods on newick tree"""

    cmd = f'TreeCluster.py -i {filename} -t {threshold} -m {method}'
    tmp = subprocess.check_output(cmd, shell=True)
    from io import StringIO
    s=str(tmp,'utf-8')
    data = StringIO(s)
    cl = pd.read_csv(data, sep='\t')
    cl.ClusterNumber+=1
    cl.ClusterNumber.replace(-1,0)
    return cl.set_index('SequenceName')

def get_treecluster_levels(tree, levels=None, method='avg_clade'):
    """Cluster a distance matrix at different threshold levels.
    Args:
        tree: newick tree file
        levels: threshold levels, a list of 1 or more levels, optional
    returns:
        a dataframe of cluster labels for each sample n the matrix
    """

    if levels == None:
        levels = [500,200,100,50,20,12,7,3]
    df = pd.DataFrame()#index=S.index)
    for t in levels:
        cl = treecluster(tree, t, method)
        #print (cl)
        df['snp'+str(t)] = cl
        cl['level'] = t
    return df

def hdbscan_cluster(distance_matrix, min_cluster_size=1, min_samples=None, alpha=1.0):
    """
    Clusters a distance matrix using HDBSCAN.

    Parameters:
        distance_matrix (pandas.DataFrame): N x N distance matrix where N is the number of samples
        min_cluster_size (int): the minimum size of clusters to be returned (default=5)
        min_samples (int or None): the number of samples in a neighborhood for a point to be considered a core point.
            If None, it will be set to min_cluster_size (default=None)
        alpha (float): the relative importance of stability vs. size in clustering (default=1.0)

    Returns:
        list: a list of cluster labels for each sample
        dict: a dictionary of the clusters that each sample is in
    """

    import hdbscan
    if min_samples is None:
        min_samples = min_cluster_size

    clusterer = hdbscan.HDBSCAN(min_cluster_size=2,metric='precomputed',alpha=alpha)
    cluster_labels = clusterer.fit_predict(distance_matrix.values)

    # create dictionary of clusters
    clusters = {}
    for i, label in enumerate(cluster_labels):
        if label != -1:
            if label in clusters:
                clusters[label].append(distance_matrix.index[i])
            else:
                clusters[label] = [distance_matrix.index[i]]

    return list(cluster_labels), clusters

def dm_cluster(distance_matrix, t, prev_clusters=None, linkage='average'):
    """
    Given a Pandas dataframe distance matrix and a distance threshold t, finds the clusters of the
    samples such that all members of a cluster are within t of each other.

    Parameters:
        distance_matrix (pandas.DataFrame): N x N distance matrix where N is the number of samples
        t (float): distance threshold
        previous_clusters (dataframe): previous cluster labels (optional)
        linkage: linkage method - 'single', 'complete', 'average'

    Returns:
        list: a list of cluster labels for each sample
        dataframe: a dataframe of the clusters that each sample is in
    """

    from sklearn.cluster import AgglomerativeClustering
    if t==0: # can't cluster at zero distance
        t=1e-10
    clustering = AgglomerativeClustering(distance_threshold=t, n_clusters=None,
                                         linkage=linkage, metric='precomputed').fit(distance_matrix)
    labels = clustering.labels_+1
    clusters = pd.DataFrame(labels,columns=['cluster'],index=distance_matrix.index)
    #print (clusters)
    if prev_clusters is not None:
        labels,clusters = reassign_clusters(clusters, labels, prev_clusters)
    return labels, clusters

def reassign_clusters(clusters, labels, prevclusters):
    """Re-assign cluster names to match previous ones."""

    new = []
    df = clusters.merge(prevclusters,left_index=True,right_index=True)
    cnt = max(prevclusters.cluster)+1
    for c,g in clusters.groupby('cluster'):
        #print (g)
        #most common cluster number from prev - in case they are mixed
        f = df.loc[df.index.isin(g.index)]
        if len(f)>0:
            #nc = f.iloc[0].cluster_y
            nc = f.cluster_y.value_counts().index[0]
            #print (f)
        else:
            nc = cnt
            cnt+=1
        #assign prev cluster instead of new one
        g['cluster'] = nc
        #print(g)
        new.append(g)

    final = pd.concat(new)
    final = final.loc[clusters.index]
    #print('curr labels:',labels)
    newlabels=list(final.cluster)
    return newlabels, final

def get_cluster_levels(S, cluster_members=None, linkage='average',
                       levels=None):
    """Cluster a distance matrix at different threshold levels.
    Args:
        S: snp distance matrix
        cluster_members (dataframe): previous sets of clusters at each level
        linkage: linkage method - 'single', 'complete', 'ward'
        levels: threshold levels, a list of 1 or more levels, optional
    returns:
        a dataframe of cluster labels for each sample n the matrix and
        a dataframe of cluster labels in long form
    """

    if levels == None:
        levels = [500,200,100,50,20,12,7,3]
    df = pd.DataFrame(index=S.index)
    clusts=[]
    for t in levels:
        if cluster_members is not None:
            prev = cluster_members[cluster_members.level==t]
        else:
            prev = None
        labels,cl = dm_cluster(S, t, prev, linkage=linkage)
        df['snp'+str(t)] = labels
        cl['level'] = t
        clusts.append(cl)
    clusts = pd.concat(clusts)
    return df, clusts

def find_reference_sample(clade, snp_distances, cl):
    """Find ref sample"""

    col='snp3'
    #print (clade)
    clade_samples = cl[cl[col] == clade].index.tolist()
    min_avg_distance = float('inf')
    reference_sample = None
    for sample in clade_samples:
        avg_distance = snp_distances.loc[clade_samples, sample].mean()
        if avg_distance < min_avg_distance:
            min_avg_distance = avg_distance
            reference_sample = sample
    return reference_sample

def generate_short_code(input_string):
    """Hexadecimal short code for input string"""

    import hashlib
    hash_object = hashlib.sha1(input_string.encode())
    hex_dig = hash_object.hexdigest()
    short_code = hex_dig[:8]
    return short_code

def generate_strain_names(cl):
    """Generate strain names from cluster levels
    Args:
        cl (dataframe): cluster level labels
    Returns:
        new strain names in a dataframe
    """

    col = 'snp3'
    col1 = 'snp7'
    col2 = 'snp12'
    col2 = 'snp50'
    col3 = 'snp200'
    col4 = 'snp500'
    # Iterate over each clade and assign a reference sample and strain names to each sample
    new = []
    for cluster,df in cl.groupby(col):
        # Find the reference sample for this clade
        #reference_sample = find_reference_sample(cluster, snpdist, cl)
        # Assign a strain name to each sample in the clade
        #df = meta[meta[col] == clade].copy()
        s=1
        df=df.replace(-1,'X')
        for id, sample in df.iterrows():
            #print (sample)
            #fs = f'{s:04d}'
            strain_name = f"ST-{sample[col4]}-{sample[col3]}-{sample[col2]}-{sample[col1]}-{sample[col]}"
            code = generate_short_code(strain_name)
            #if id == reference_sample:
            #    strain_name += "-ref"
            s+=1
            df.loc[id,'strain_name'] = strain_name
            df.loc[id,'code'] = code
        #cladecount+=1
        #print (df[['strain_name','code']])
        new.append(df)
    new = pd.concat(new)
    return new

def nonredundant_samples(df, col='snp3'):
    """Get non redundant samples at a specific threshold from cluster table"""

    df['dup'] = df[col].replace(-1,np.nan)
    x = df[(~df.duplicated('dup')) | (df.dup.isnull())]
    return x