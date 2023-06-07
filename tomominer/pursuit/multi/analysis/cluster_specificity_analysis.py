#!/usr/bin/env python



'''
use this to inspect alignment scores etc to see why cluster_removal_according_to_center_matching_specificity() sometimes cannot mark some redundant clusters as none-specific
'''




# collect alignment scores of alignment of subtomgrams of one cluster (self_clus) to its own average (self_tem) and to another average (other_tem), list these scores then perform wilcox test
if __name__ == '__main__':

    import sys
    data_dir = sys.argv[1]
    self_clus = int(sys.argv[2])

    import os
    import pickle
    with open(os.path.join(data_dir, 'cluster_average_select.pickle')) as f: cas = pickle.load(f)
    with open(os.path.join(data_dir, 'align_template.pickle')) as f: al = pickle.load(f)


    self_tem = cas['selected_templates'][self_clus]['subtomogram'];     print 'this cluster average', self_tem

    for other_clus in cas['selected_templates']:
        if other_clus == self_clus:     continue

        print other_clus, '\t\t\t\t',

        other_tem = cas['selected_templates'][other_clus]['subtomogram'];        print other_tem
 
        self_clus_member = set(_['subtomogram'] for _ in cas['tk_info'][self_tem]['data_json'])


        self_scores = []
        other_scores = []
        for r in al:
            r = r.result
            if r['vol_key']['subtomogram'] not in self_clus_member:    continue

            self_scores.append(r['align'][self_clus]['score'])
            other_scores.append(r['align'][other_clus]['score'])


        #for i in range(len(self_scores)):       print '%0.3f'%(self_scores[i] - other_scores[i]), 

        import numpy as N
        import scipy.stats as SCS
        print SCS.wilcoxon(other_scores, self_scores), N.median(self_scores) - N.median(other_scores)


