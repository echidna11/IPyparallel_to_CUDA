#!/usr/bin/env python



'''
use this to inspect alignment scores etc to see why cluster_removal_according_to_center_matching_specificity() sometimes cannot mark some redundant clusters as none-specific
the alignment scores are obtained from tomominer.template alignment calculated here, not stored ones. So that we can test different alignment parameters's effect to the specificity
'''




# collect alignment scores of alignment of subtomgrams of one cluster (self_clus) to its own average (self_tem) and to another average (other_tem), list these scores then perform wilcox test
if __name__ == '__main__':

    import sys
    data_dir = sys.argv[1]
    self_clus = int(sys.argv[2])
    other_clus = int(sys.argv[3])

    import os
    import pickle
    with open(os.path.join(data_dir, 'cluster_average_select.pickle')) as f: cas = pickle.load(f)


    self_tem = cas['selected_templates'][self_clus]['subtomogram'];     print self_tem
    other_tem = cas['selected_templates'][other_clus]['subtomogram'];        print other_tem

    import tomominer.common.obj as CO
    self = CO.Object()
    
    from tomominer.io.cache import Cache
    self.cache = Cache(tmp_dir=os.getenv('TMP_DIR'))

    from tomominer.parallel.queue_master import QueueMaster
    self.runner = QueueMaster(host='localhost', port=5011)

    if True:
        stk = cas['selected_templates']
    else:
        stk = cas['selected_templates_common_frame']
 
    tk = {}
    tk[self_clus] = stk[self_clus]
    tk[other_clus] = stk[other_clus]

    tasks = []
    for r in cas['tk_info'][self_tem]['data_json']:   
        tasks.append(self.runner.task( module='tomominer.pursuit.multi.util', method='align_to_templates', kwargs={'rec':r, 'tem_keys':tk, 'align_op':{'L':36, 'with_missing_wedge':True}, 'multiprocessing':False} ))

    self_scores = []
    other_scores = []
    for r in self.runner.run__except(tasks):
        r = r.result

        self_scores.append(r['align'][self_clus]['score'])
        other_scores.append(r['align'][other_clus]['score'])

    for i in range(len(self_scores)):       print '%0.3f'%(self_scores[i] - other_scores[i]), 

    import numpy as N
    import scipy.stats as SCS
    print SCS.wilcoxon(other_scores, self_scores), N.median(self_scores) - N.median(other_scores)


