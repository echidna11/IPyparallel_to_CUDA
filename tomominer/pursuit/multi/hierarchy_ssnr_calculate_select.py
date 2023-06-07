

'''
after dimension reduction
first perform hiearachical clustering, and record the member indicis to JSON list
then for every cluster node that contains more than XX members, assign a ssnr=None field in order to save SSNR statistics
then repeatly start from tomominer.clusters whose children either have SSNR assigned or do not have ssnr field, calculate its statistics of SSNR, calculate its SSNR and FSC and assign to the cluster


the calculation can be done in parallel, either through computer cluster or through multiprocessing.
* for computer cluster version, before assign uuid field inside ssnr field, in order to store and load SSNR statistics. After finishing the whole calculation, delete all SSNR statistics files according to uuids
* For the multi-processing version, the ssnr statistics can be directly returned by the called function, instead of writing to disk.

'''

import os
import sys


import scipy.spatial.distance
import scipy.cluster.hierarchy
import tomominer.cluster.hierarchy.hierarchy as CH
import fastcluster
# perform hierachical clustering, get tree representation t, and collect information i
def hierarchical_clustering(x):

    c = fastcluster.average(scipy.spatial.distance.pdist(x))
    t = scipy.cluster.hierarchy.to_tree(c)
    i = CH.collect_cluster_info(tree=t, x=x)

    return {'info':i, 'root_id':t.id, 'clusters':c}




from collections import deque
'''
breadth first traverse from root node
parameters:     hi: return of CH.collect_cluster_info(),    root: root id
'''

def hierarchy_traverse(hi, root):
    q = deque()            # node queue
    q.append(root)

    while len(q) > 0:
        n = q.popleft()
        yield n

        if 'left' in hi[n]:      q.append(hi[n]['left'])
        if 'right' in hi[n]:     q.append(hi[n]['right'])    


import pickle
# get ssnr statistics and delete no need fields to save memory
def ssnr_stat__load(sd):

    if 'key' in sd:
        key = sd['key']             ;               assert      key is not None
        with open(key, 'rb') as f:    ssnr = pickle.load(f)

    else:

        ssnr = {}
        ssnr['sum'] = sd['sum'].copy()
        ssnr['prod_sum'] = sd['prod_sum'].copy()
        ssnr['mask_sum'] = sd['mask_sum'].copy()

        del sd['sum']
        del sd['prod_sum']
        del sd['mask_sum']


    return ssnr



import copy
# add up ssnrs and return the sum
def ssnr_statistics_combine(s1, s2):
    s = None

    if s1 is not None:
        s = copy.deepcopy(s1)

    if s2 is not None:
        if s is None:
            s = copy.deepcopy(s2)
        else:
            s['sum'] += s2['sum']
            s['prod_sum'] += s2['prod_sum']
            s['mask_sum'] += s2['mask_sum']

    return s


import tomominer.statistics.ssnr as SS
def ssnr_stat_local(self, dj, return_key, segmentation_tg_op=None):
    s = SS.var__local(self=self, data_json=dj, return_key=return_key, segmentation_tg_op=segmentation_tg_op)

    s['sum'] = s['sum'][0]
    s['prod_sum'] = s['prod_sum'][0]
    s['mask_sum'] = s['mask_sum'][0]

    return s



# calculate ssnr and its statistics for a single node, delete the statistics of its children to save storage
# parameters:       n_id: node id         hi: node information        dj: the dictionary of data_json that only contains the data of nl_id and nr_id
def hierarchy_ssnr_calculate__single_node(self, n_id, hi, djd, return_key=False, segmentation_tg_op=None, ssnr_op=None):


    if 'left' in hi[n_id]:
        nl_id = hi[n_id]['left']
        if 'ssnr' in hi[nl_id]:
            assert      hi[nl_id]['ssnr'] is not None
            ssnr_l = ssnr_stat__load(hi[nl_id]['ssnr'])
        else:
            ssnr_l = ssnr_stat_local(self, djd[nl_id], return_key=False, segmentation_tg_op=segmentation_tg_op)           # if the child does not have ssnr statistics, calculate them
            

    if 'right' in hi[n_id]:
        nr_id = hi[n_id]['right']
        if 'ssnr' in hi[nr_id]:
            assert      hi[nr_id]['ssnr'] is not None
            ssnr_r = ssnr_stat__load(hi[nr_id]['ssnr'])
        else:
            ssnr_r = ssnr_stat_local(self, djd[nr_id], return_key=False, segmentation_tg_op=segmentation_tg_op)           # if the child does not have ssnr statistics, calculate them


    ssnr_stat = None
    ssnr_stat = ssnr_statistics_combine(ssnr_stat, ssnr_l)
    ssnr_stat = ssnr_statistics_combine(ssnr_stat, ssnr_r)


    ssnr = SS.ssnr__given_stat(sum_v=ssnr_stat['sum'], prod_sum=ssnr_stat['prod_sum'], mask_sum=ssnr_stat['mask_sum'], op=ssnr_op)



    if return_key:

        re_key = self.cache.save_tmp_data(ssnr_stat, fn_id=self.task.task_id)
        assert      re_key is not None              # if two workers A and B are working on the same task, and the result has been already saved by worker B, then when A calls self.cache.save_tmp_data() will return None. In this case, the assertion will result an exception so that hierarchy_ssnr_calculate__single_node() will fail and worker A will go ahead for other pending tasks
        ssnr['key'] = re_key

    else:

        ssnr['sum'] = ssnr_stat['sum']
        ssnr['prod_sum'] = ssnr_stat['prod_sum']
        ssnr['mask_sum'] = ssnr_stat['mask_sum']

    return {'id':n_id, 'ssnr':ssnr}



def collect__info_and_data(n_id, hi, dj):
    
    hi_t = {}
    hi_t[n_id] = hi[n_id]

    djd = {}

    if 'left' in hi[n_id] is not None:      
        idt = hi[n_id]['left']
        hi_t[idt] = hi[idt]

        if 'ssnr' not in hi[idt]:
            # this means this child need to calculate ssnr statistics using data_json
            djd[idt] = [dj[_] for _ in hi[idt]['nodes']]

    del idt

    if 'right' in hi[n_id]:      
        idt = hi[n_id]['right']
        hi_t[idt] = hi[idt]

        if 'ssnr' not in hi[idt]:
            djd[idt] = [dj[_] for _ in hi[idt]['nodes']]


    del idt

    return {'hi':hi_t, 'djd':djd}



# clean up ssnr statistics
def clean_up_ssnr_stat(hi_t):
    if 'ssnr' not in hi_t:          return False
    if hi_t['ssnr'] is None:        return False

    do_clean = False

    if 'sum' in hi_t['ssnr']:      
        del hi_t['ssnr']['sum']
        do_clean = True

    if 'prod_sum' in hi_t['ssnr']:      
        del hi_t['ssnr']['prod_sum']
        do_clean = True

    if 'mask_sum' in hi_t['ssnr']:      
        del hi_t['ssnr']['mask_sum']
        do_clean = True

    return do_clean



import tomominer.common.obj as CO

# use multiprocessing for parallel
# parameters:       n: node         hi: node information        dj: data_json
def hierarchy_ssnr_calculate__multiprocessing(self, ns, hi, dj, segmentation_tg_op=None):


    clear_count = 0
    for i in hi:
        if 'parent' not in hi[i]:  continue
        pi = hi[i]['parent']
        if 'ssnr' not in hi[pi]:    continue
        if hi[pi]['ssnr'] is None:   continue
        
        # if this node's parent already has ssnr statistics calculated, remove this node's ssnr statistics in order to save memory
        if clean_up_ssnr_stat(hi[i]):        clear_count += 1

    print '\r\t\t\t\t\t\t\t\t', clear_count, 'ssnr statistics cleaned up',

    self_trim = CO.Object();       self_trim.cache = self.cache        # runner and pool cannot be passed through multiprocessing!!!
    
    pool_re = []
    for n_id in ns:
        cid = collect__info_and_data(n_id, hi, dj)
        #hierarchy_ssnr_calculate__single_node(self=self_trim, n_id=n_id, hi=cid['hi'], djd=cid['djd'], return_key=False)
        pool_re.append(    self.pool.apply_async(   func=hierarchy_ssnr_calculate__single_node, kwds={'self':self_trim, 'n_id':n_id, 'hi':cid['hi'], 'djd':cid['djd'], 'return_key':False, 'segmentation_tg_op':segmentation_tg_op}    )      )
        

    for r in pool_re:
        r = r.get(999999)
        hi[r['id']]['ssnr'] = r['ssnr']



        

def hierarchy_ssnr_calculate__computer_cluster(self, ns, hi, dj, segmentation_tg_op=None):
    tasks = []
    for n_id in ns:
        cid = collect__info_and_data(n_id, hi, dj)
        tasks.append(    self.runner.task(module='tomominer.pursuit.multi.hierarchy_ssnr_calculate_select', method='hierarchy_ssnr_calculate__single_node', kwargs={'n_id':n_id, 'hi':cid['hi'], 'djd':cid['djd'], 'return_key':True, 'segmentation_tg_op':segmentation_tg_op})    )

    for r in self.runner.run__except(tasks):
        r = r.result
        hi[r['id']]['ssnr'] = r['ssnr']



import multiprocessing
from multiprocessing.pool import Pool as Pool
import gc
# parallel calculate ssnr for each cluster
# parameters:   hc: hieratchical clustering information,        dj: data_json
def hierarchy_ssnr_calculate(self, hc, dj, cluster_size_min=0, cluster_size_max=float('inf'), multiprocessing_factor=1.0, segmentation_tg_op=None):

    print 'hierarchy_ssnr_calculate', cluster_size_max

    if self.pool is None:        self.pool = Pool()

    hi = hc['info']


    # assign empty ssnr field for all large clusters
    for n_id in hi:
        if (len( hi[n_id]['nodes'] ) >= cluster_size_min) and (len( hi[n_id]['nodes'] ) <= cluster_size_max):
            hi[n_id]['ssnr'] = None
        else:
            assert      'ssnr' not in hi[n_id]



    process_pass = 0
    while True:

        # collect cluster nodes that need to calculate ssnr
        ns = []
        for n_id in hi:
            if 'ssnr' not in hi[n_id]:      continue        # this cluster is not eligible to calculate ssnr
            if hi[n_id]['ssnr'] is not None:        continue            # this cluster already has its ssnr calculated

            if ('left' in hi[n_id]) and ('ssnr' in hi[hi[n_id]['left']]) and (hi[hi[n_id]['left']]['ssnr'] is None):      continue        # this node's child need to calculate SSNR statistics
            if ('right' in hi[n_id]) and ('ssnr' in hi[hi[n_id]['right']]) and (hi[hi[n_id]['right']]['ssnr'] is None):      continue        # this node's child need to calculate SSNR statistics

            ns.append(n_id)

        print '\r', 'pass', process_pass, 'processing', len(ns), 'clusters        ',
        sys.stdout.flush()

        if len(ns) == 0:    break       # no more to process

        if len(ns) == 1:
            # in this case, directly call the ssnr calculation rotine for the single node
            cid = collect__info_and_data(ns[0], hi, dj)
            hi[ns[0]]['ssnr'] = hierarchy_ssnr_calculate__single_node(self=self, n_id=ns[0], hi=cid['hi'], djd=cid['djd'], return_key=False, segmentation_tg_op=segmentation_tg_op)['ssnr']

        elif len(ns) <= (multiprocessing.cpu_count() * multiprocessing_factor):
            # if there are no many eligible nodes, use multiprocessing
            hierarchy_ssnr_calculate__multiprocessing(self=self, ns=ns, hi=hi, dj=dj, segmentation_tg_op=segmentation_tg_op)

        else:
            # if there are so many clusters, use the computer cluster for calculation
            hierarchy_ssnr_calculate__computer_cluster(self=self, ns=ns, hi=hi, dj=dj, segmentation_tg_op=segmentation_tg_op)

        process_pass += 1



    # clean up ssnr statistics
    for i in hi:        clean_up_ssnr_stat(hi[i])

    # clean up temporary files
    for i in hi:
        if 'ssnr' not in hi[i]: continue
        if 'key' not in hi[i]['ssnr']: continue

        try:
            os.remove(hi[i]['ssnr']['key'])
        except:
            pass

        del hi[i]['ssnr']['key']

    for i in hi:
        if 'ssnr' not in hi[i]:     continue
        assert      hi[i]['ssnr'] is not None


    self.pool.close()
    self.pool = None


# order all clusters according to decrease of ssnr, and select the largest set of disjoint clusters with highest ssnr
def hierarchy_ssnr_cluster_select(hi):
    ids = sorted([_ for _ in hi if 'ssnr' in hi[_]], key=lambda _:(-hi[_]['ssnr']['fsc'].sum()) )     # sort cluster ids according to decrease of FSC

    ids_selected = []
    nodes = set()
    for i in ids:
        if len(nodes.intersection(hi[i]['nodes'])) > 0:   continue
        ids_selected.append(i)
        nodes.update(hi[i]['nodes'])

    return ids_selected


# convert selected clusters to labels that is consistant with data_json
def cluster_select_label(hi, ids, dj_len):
    labels = [-1] * dj_len

    current_label = 0
    id_to_label = {}
    label_to_id = {}

    for id_t in ids:
        for i in hi[id_t]['nodes']:
            assert labels[i] < 0
            labels[i] = current_label
        
        id_to_label[id_t] = current_label
        label_to_id[current_label] = id_t

        current_label += 1

    return {'labels':labels, 'id_to_label':id_to_label, 'label_to_id':label_to_id}






