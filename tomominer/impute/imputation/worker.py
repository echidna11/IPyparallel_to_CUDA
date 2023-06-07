from worker import ClassifyWorker

class ClassifyWithImputationWorker(ClassifyWorker) :
    
    def __init__(self, host, port):
        ClassifyWorker.__init__(self, host, port) #, level=logging.WARN)

    
    # todo: make a distributed k-means, where subtomograms are loaded, then imputed
    def cluster(self,mat_structure, dims, k):

        from scipy.cluster.vq import kmeans, vq

        mat = np.vstack( np.hstack(np.load(_) for _ in row) for row in mat_structure )


        centroids, _ = kmeans(red, k)
        labels,    _ = vq(red, centroids)

        return labels

    # performing k-means clustering at local node, for each subtomogram, find its closest cluster center, then assign the cluster label, and calculate local sum
    def cluster_local(self, vmal, cc) :
        for vk, mk, ang, loc in vmal:
            vol = self.get_mrc(vk)



    # mxu: todo: add interpolation, and remove mask 
    def align_to_templates(self, vmal, vm_tem, L):

        v1_key, m1_key, ang, loc = vmal

        v1 = self.get_mrc(v1_key)
        m1 = self.get_mrc(m1_key)

        best = type('', (), {})()       # create an object for assigning dynamic attributes
        best.subtomo_key = v1_key
        best.subtomo_msk_key = m1_key

        best.score = float('-inf')

        for i,(v2_key,m2_key) in enumerate(vm_tem):
            v2 = self.get_mrc(v2_key)
            m2 = self.get_mrc(m2_key)

            ea = np.array([0.0,0.0,0.0])
            dx = np.array([0.0,0.0,0.0])

            score = tomo.combined_search_py(v2, m2, v1, m1, L, dx, ea)

            if score > best.score:
                best.score = score
                best.ang = ea
                best.loc = dx
                best.tem_key = v2_key
                best.tem_msk_key = m2_key

        
        return best



