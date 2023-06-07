
import numpy as np

def kmeans(vols, k):
    import scipy.cluster.vq as vq
    
    vol_mat = np.vstack( np.transpose(vol.flatten())  for vol in vols )  # first convert vols to row vectors, and forma data matrix

    print vol_mat.shape

    centroids, _ = vq.kmeans(vol_mat, k)
    labels,    _ = vq.vq(vol_mat, centroids)

    return labels

def cluster_average(vols, labels):
   
   vol_avg = {}

   for lbl in set(labels):
       vol_avg[lbl] = sum(vol[labels == lbl, :, :]) / sum(labels == lbl) 
    
   return vol_avg 

