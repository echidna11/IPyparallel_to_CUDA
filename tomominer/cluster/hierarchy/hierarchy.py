

# functions for hierarchical clustering
import copy
import numpy as N



# collect all clusters inside the hierarchy, for each cluster, record its tightness and nodes
# input, hierarchical clustering encoded as a linkage matrix, outputed from scipy.cluster.hierarchy.linkage()
def collect_cluster_info(tree, x=None):
    info_d = {}

    for n in hierarchy_traverse(tree):
        info_d[n.id] = {}

    for n in hierarchy_traverse(tree):
        if n.left is not None:
            info_d[n.id]['left'] = n.left.id
            info_d[n.left.id]['parent'] = n.id

        if n.right is not None:      
            info_d[n.id]['right'] = n.right.id
            info_d[n.right.id]['parent'] = n.id


    for n in hierarchy_traverse(tree):
        if not n.is_leaf():
            info = {'size':n.get_count(), 'dist':n.dist, 'nodes':[_.id for _ in hierarchy_traverse(n) if _.is_leaf()]}        # n.get_dist() outputs cluster tightness

            if x is not None:                info['centroid'] = x[info['nodes'], :].mean(axis=0)
        else:
            info = {'size':1, 'dist':n.dist, 'nodes':[n.id]}        # n.get_dist() outputs cluster tightness

            if x is not None:                info['centroid'] = x[n.id, :]

        info_d[n.id] = dict(info_d[n.id].items() + info.items())
    
    return info_d



    
from collections import deque
# breadth first traverse from root node
def hierarchy_traverse(root):
    q = deque()            # node queue
    q.append(root)

    while len(q) > 0:
        n = q.popleft()
        yield n

        if n.left is not None:      q.append(n.left)
        if n.right is not None:     q.append(n.right)    




def cluster_filter(info, min_size=None):
    info_f = {}
    for id_t in info:
        if min_size is not None:
            if info[id_t]['size'] < min_size:   continue
        
        info_f[id_t] = info[id_t]
        
    
    return info_f            



#----------------------------------------------
# calculate statistics

'''
peform width first search to obtain the lists of cluster root ids, given different cutoffs, just like we cut on a hierarchical clustering tree
input, r: root id of hierarchical clustering encoded by output of collect_cluster_info(),     hi: tree encoded by collect_cluster_info()
'''
def cluster_number_vs_dist_cutoff(r, hi):
    level_records = []
    
    from Queue import PriorityQueue

    q = PriorityQueue()
    q.put((-hi[r]['dist'], r))
    level_records.append({'dist_cutoff':hi[r]['dist'], 'cluster_num':q.qsize()})

    while not q.empty():
        n = q.get()[1]
        if hi[n]['dist'] <= 0:      continue            # do not process sigletons

        if 'left' in hi[n]:
            nt = hi[n]['left']
            q.put((-hi[nt]['dist'], nt))

        if 'right' in hi[n]:
            nt = hi[n]['right']
            q.put((-hi[nt]['dist'], nt))

        level_records.append({'dist_cutoff':hi[n]['dist'], 'cluster_num':q.qsize()})
        #print level_records[-1]['dist_cutoff'], level_records[-1]['cluster_num'], '\t',
    
    return level_records


'''
given a distance cutoff, get corresponding clusters in the hierarchy
'''
def cluster_cut(r, hi, dist_cutoff):
    q = deque()            # node queue
    q.append(r)

    sel_clus = []
    while len(q) > 0:
        n = q.popleft()

        if hi[n]['dist'] <= dist_cutoff:
            sel_clus.append(n)
            continue            # in such case, do not further expand / search subclusters

        if 'left' in hi[n]:     q.append(hi[n]['left'])
        if 'right' in hi[n]:     q.append(hi[n]['right'])

    return sel_clus



#---------------------------------------
# functions for getting optimal clusters

# perform hierarchical clustering given a distance matrix (in square form)
def cluster_dist_mat(d):
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import squareform
    return linkage(squareform(d))



# calculate Silhouette for each level and see which is optimal
# d is the distance matrix in square form, z is the return from scipy.cluster.hierarchy.linkage()
def optimal_silhouette(d, z, dist_threhold_max= float('inf')):

    from scipy.cluster.hierarchy import fcluster
    from sklearn.metrics import silhouette_score

    best = {'score': (-float('inf'))}
    for i in range(z.shape[0] - 1):
        t = z[i, 2]
        if t > dist_threhold_max:       continue

        lbl = fcluster(z, t, criterion='distance')
        ss = silhouette_score(d, lbl, metric='precomputed')

        if ss > best['score']:
            best['score'] = ss
            best['labels'] = lbl
            best['threshold'] = t

    return best

'''
given a cluster labelling l0, get the cutoff that is most consistent with l0
parameters:     z is the return from scipy.cluster.hierarchy.linkage(),         find_min_dist=True: if there are multiple best choices, find the one with minimal cutoff,   find_min_dist=False: find one with maximal cutoff.          consistency_matric={'ari', 'mi'}: ari:Adjusted Rand index,      mi: Mutual Information,     ami: Adjusted Mutual Information,       nmi: Normalized Mutual Information
'''
def optimal_consistent_cut(z, l0, find_min_dist=True, consistency_matric='nmi', verbose=False):
    l0 = N.array(l0)
    ind = N.isfinite(l0)

    dist_s = z[:, 2].flatten()
    dist_s = dist_s[dist_s > 0]
    dist_s = dist_s.tolist()

    if find_min_dist:
        dist_s = sorted(dist_s, reverse=True)
    else:
        dist_s = sorted(dist_s)

    from scipy.cluster.hierarchy import fcluster
    from sklearn import metrics
    best = None
    for t in dist_s:
        l = fcluster(z, t, criterion='distance')

        if consistency_matric == 'ari':
            c = metrics.adjusted_rand_score(l0[ind], l[ind])
        elif consistency_matric == 'mi':
            c = metrics.mutual_info_score(l0[ind], l[ind])
        elif consistency_matric == 'ami':
            c = metrics.adjusted_mutual_info_score(l0[ind], l[ind])
        elif consistency_matric == 'nmi':
            c = metrics.normalized_mutual_info_score(l0[ind], l[ind])
        else:
            raise   Exception('consistency_matric')
        
        if (best is None) or (best['c'] <= c):
            best = {'c':c, 't':t}
            if verbose:   print best

    return best




'''
# testing commands


import os,sys
sys.path.append(os.path.join( os.getenv('HOME'), 'ln/tomo/py' ))

import numpy as N
x = N.random.random( (10, 3) )

from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
d = squareform(pdist(x))

import tomominer.cluster.hierarchy as H
z = H.cluster_dist_mat(d)

bs = H.optimal_silhouette(d, z)
print bs



from scipy.cluster.hierarchy import to_tree
zt = to_tree(z)
ci = H.collect_cluster_info(zt)



'''


