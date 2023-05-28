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

snp200_cmap = {0:'darkgreen',1:'coral',2:'dodgerblue',3:'crimson',4:'mediumpurple',
                5:'orange',6:'pink',7:'gold',8:'cyan',9:'beige',10:'red',
                11:'brown',12:'gray',13:'blue',14:'lightgreen',15:'Aquamarine',
                16:'DarkCyan',17:'IndianRed'}

def prep(tree, support, resolve_polytomies=True, suppress_unifurcations=True):

    if resolve_polytomies:
        tree.resolve_polytomies()
    if suppress_unifurcations:
        tree.suppress_unifurcations()
    leaves = set()
    for node in tree.traverse_postorder():
        if node.edge_length is None:
            node.edge_length = 0
        node.DELETED = False
        if node.is_leaf():
            leaves.add(str(node))
        else:
            try:
                node.confidence = float(str(node))
            except:
                node.confidence = 100. # give edges without support values support 100
            if node.confidence < support: # don't allow low-support edges
                node.edge_length = float('inf')
    return leaves

def cut(node):
    """cut out the current node's subtree (by setting all nodes' DELETED to True) 
    and return list of leaves"""

    from queue import PriorityQueue,Queue
    cluster = list()
    descendants = Queue(); descendants.put(node)
    while not descendants.empty():
        descendant = descendants.get()
        if descendant.DELETED:
            continue
        descendant.DELETED = True
        descendant.left_dist = 0; descendant.right_dist = 0; descendant.edge_length = 0
        if descendant.is_leaf():
            cluster.append(str(descendant))
        else:
            for c in descendant.children:
                descendants.put(c)
    return cluster

def single_linkage_cut(tree,threshold,support):
    """single-linkage clustering using Metin's cut algorithm"""

    leaves = prep(tree,support)
    clusters = list()

	# find closest leaf below (dist,leaf)
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.min_below = (0,node.label)
        else:
            node.min_below = min((c.min_below[0]+c.edge_length,c.min_below[1]) for c in node.children)

    # find closest leaf above (dist,leaf)
    for node in tree.traverse_preorder():
        node.min_above = (float('inf'),None)
        if node.is_root():
            continue
        # min distance through sibling
        for c in node.parent.children:
            if c != node:
                dist = node.edge_length + c.edge_length + c.min_below[0]
                if dist < node.min_above[0]:
                    node.min_above = (dist,c.min_below[1])
        # min distance through grandparent
        if not c.parent.is_root():
            dist = node.edge_length + node.parent.min_above[0]
            if dist < node.min_above[0]:
                node.min_above = (dist,node.parent.min_above[1])

    # find clusters
    for node in tree.traverse_postorder(leaves=False):
        # assume binary tree here (prep function guarantees this)
        l_child,r_child = node.children
        l_dist = l_child.min_below[0] + l_child.edge_length
        r_dist = r_child.min_below[0] + r_child.edge_length
        a_dist = node.min_above[0]
        bad = [0,0,0] # left, right, up
        if l_dist + r_dist > threshold:
            bad[0] += 1; bad[1] += 1
        if l_dist + a_dist > threshold:
            bad[0] += 1; bad[2] += 1
        if r_dist + a_dist > threshold:
            bad[1] += 1; bad[2] += 1
        # cut either (or both) children
        for i in [0,1]:
            if bad[i] == 2:
                cluster = cut(node.children[i])
                if len(cluster) != 0:
                    clusters.append(cluster)
                    for leaf in cluster:
                        leaves.remove(leaf)
        # cut above (equals cutting me)
        if bad[2] == 2: # if cutting above, just cut me
            cluster = cut(node)
            if len(cluster) != 0:
                clusters.append(cluster)
                for leaf in cluster:
                    leaves.remove(leaf)
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters

# min_clusters_threshold_max, but all clusters must define a clade
def min_clusters_threshold_max_clade(tree,threshold,support):
    leaves = prep(tree, support, resolve_polytomies=False)

    # compute leaf distances and max pairwise distances
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.leaf_dist = 0; node.max_pair_dist = 0
        else:
            node.leaf_dist = float('-inf'); second_max_leaf_dist = float('-inf')
            for c in node.children: # at least 2 children because of suppressing unifurcations
                curr_dist = c.leaf_dist + c.edge_length
                if curr_dist > node.leaf_dist:
                    second_max_leaf_dist = node.leaf_dist; node.leaf_dist = curr_dist
                elif curr_dist > second_max_leaf_dist:
                    second_max_leaf_dist = curr_dist
            node.max_pair_dist = max([c.max_pair_dist for c in node.children] + [node.leaf_dist + second_max_leaf_dist])

    # perform clustering
    q = Queue(); q.put(tree.root); roots = list()
    while not q.empty():
        node = q.get()
        if node.max_pair_dist <= threshold:
            roots.append(node)
        else:
            for c in node.children:
                q.put(c)

    # if verbose, print the clades defined by each cluster
    #if VERBOSE:
    #    for root in roots:
    #        print("%s;" % root.newick(), file=stderr)
    return [[str(l) for l in root.traverse_leaves()] for root in roots]

def get_treecluster_labels(clusters, label='label'):

    labels = []
    cluster_num = 1
    for cluster in clusters:
        if len(cluster) == 1:
            l = list(cluster)[0]
            labels.append([l,-1])
        else:
            for l in cluster:
                labels.append([l,cluster_num])
            cluster_num += 1

    cl = pd.DataFrame(labels,columns=['id',label])
    return cl

def treecluster(filename):
    """Run treecluster methods on newick tree"""

    from treeswift import read_tree_newick

    infile = open(filename)
    for line in infile:
            if isinstance(line,bytes):
                l = line.decode().strip()
            else:
                l = line.strip()

    tree = read_tree_newick(l)
    clusters = min_clusters_threshold_max_clade(tree,18,0)
    cl = get_treecluster_labels(clusters)
    return cl

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

def dm_cluster(distance_matrix, t, prev_clusters=None):
    """
    Given a Pandas dataframe distance matrix and a distance threshold t, finds the clusters of the
    samples such that all members of a cluster are within t of each other.
    
    Parameters:
        distance_matrix (pandas.DataFrame): N x N distance matrix where N is the number of samples
        t (float): distance threshold
        previous_clusters (dataframe): previous cluster labels (optional)
        
    Returns:
        list: a list of cluster labels for each sample
        dataframe: a dataframe of the clusters that each sample is in
    """
    
    from sklearn.cluster import AgglomerativeClustering
    clustering = AgglomerativeClustering(distance_threshold=t, n_clusters=None, 
                                         linkage='complete', metric='precomputed').fit(distance_matrix)
    labels = clustering.labels_+1
    clusters = pd.DataFrame(labels,columns=['cluster'],index=distance_matrix.index)
    if prev_clusters is not None:
        labels,clusters = reassign_clusters(clusters, labels, prev_clusters)
    return labels, clusters

def reassign_clusters(clusters, labels, prevclusters):
    """Re-assign cluster names to match previous ones."""

    new=[]
    df=clusters.merge(prevclusters,left_index=True,right_index=True)
    cnt = max(prevclusters.cluster)+1
    for c,g in clusters.groupby('cluster'):
        #print (c)       
        f = df.loc[df.index.isin(g.index)]
        #print (f)
        if len(f)>0:
            nc = f.iloc[0].cluster_y            
        else:
            nc = cnt
            cnt+=1
        g['cluster'] = nc
        #print(g)
        new.append(g)
        
    final = pd.concat(new)
    final = final.loc[clusters.index]
    #print('curr labels:',labels)
    newlabels=list(final.cluster)
    return newlabels, final

def get_cluster_levels(S, cluster_members=None):
    """Clusters at different thresholds.
    S: snp distance matrix  
    cluster_members (dataframe): previous sets of clusters at each level    
    """

    levels=[500,200,50,12,3]
    df=pd.DataFrame(index=S.index)
    clusts=[]
    for t in levels:
        if cluster_members is not None:
            prev = cluster_members[cluster_members.level==t]
        else:
            prev = None
        labels,cl = dm_cluster(S, t, prev)
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
            strain_name = f"ST-{sample[col4]}-{sample[col3]}-{sample[col2]}-{sample[col]}"
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