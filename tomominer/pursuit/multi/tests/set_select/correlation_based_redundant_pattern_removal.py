
'''
see
https://docs.google.com/document/d/1oXEGGmHOr_fUEhyoWeCavKfpdBjyoegBAVdMJDPdR38/edit#
'''

import os, json, copy
import cPickle as pickle
import numpy as N

def main0():
    with open(file_stat_file) as f:     fs = json.load(f)
    fs['passes'] = {int(_):fs['passes'][_] for _ in fs['passes']}

    # collect specificity, and correlation information for a particular pass
    fst = fs['passes'][pass_i]

    with open(os.path.join(fst['pass_dir'], 'cluster_average_align_common_frame__multi_pair.pickle'), 'rb') as f:   caacfmp = pickle.load(f)
    with open(os.path.join(fst['pass_dir'], 'cluster_removal_according_to_center_matching_specificity.pickle'), 'rb') as f:     cratcms = pickle.load(f)

    tka = caacfmp['tka']        # aligned selected templates
    pa = caacfmp['pa']
    wst = cratcms['wilcoxion_stat']

    for c0 in pa:
        if c0 not in wst:   continue
        for c1 in pa[c0]:
            if c1 not in wst[c0]:       continue
            print c0, c1, pa[c0][c1]['score'], wst[c0][c1]



def main():
    with open(file_stat_file) as f:     fs = json.load(f)
    fs['passes'] = {int(_):fs['passes'][_] for _ in fs['passes']}

    # collect specificity, and correlation information for a particular pass
    fst = fs['passes'][pass_i]

    cluster_info = {}
    for i in fs['passes']:
        with open(fs['passes'][i]['cluster_info_file'], 'rb') as f:     cluster_info[i] = pickle.load(f)

    with open(fs['passes'][i]['cluster_info_stat_file'], 'rb') as f:     cluster_info_stat = pickle.load(f)
    with open(fs['passes'][i]['align_template_file'], 'rb') as f:     at_ress = pickle.load(f)
    with open(fs['passes'][i]['cluster_average_select_file'], 'rb') as f:     cas_re = pickle.load(f)

    st = align_compare_stat(ci=cluster_info, cis=cluster_info_stat, al=[_.result for _ in at_ress], tk=cas_re['selected_templates'])

    with open(os.path.join(fst['pass_dir'], 'cluster_average_align_common_frame__multi_pair.pickle'), 'rb') as f:   caacfmp = pickle.load(f)
    pa = caacfmp['pa']

    # verification of consistancy of the cluster id
    for c in caacfmp['tka']:
        assert caacfmp['tka'][c]['pass_i'] == cas_re['selected_templates'][c]['pass_i']
        assert caacfmp['tka'][c]['cluster'] == cas_re['selected_templates'][c]['cluster']
        assert caacfmp['tka'][c]['id'] == c
        assert cas_re['selected_templates'][c]['id'] == c

    for c0 in pa:
        if c0 not in st:   continue
        for c1 in pa[c0]:
            #if c1 <= c0:        continue        # we assume the clusters are ordered according to FSC scores
            if c1 not in st[c0]:       continue
            print c0, c1, pa[c0][c1]['score'], (st[c0][c1]['median0'] > st[c0][c1]['median1']), st[c0][c1]['p']


    import pdb;     pdb.set_trace()



'''
this function is modified from 
util.cluster_removal_according_to_center_matching_specificity()
'''
import scipy.stats as SCS
def align_compare_stat(ci, cis, al, tk, test_sample_num_min=10):
    print 'align_compare_stat()'

    tk = copy.deepcopy(tk)
    tkd = {tk[_]['subtomogram']:_ for _ in tk}      # mapping from tomominer.template key to cluster id

    from collections import defaultdict
    st = defaultdict(dict)
    for pass_i in ci:
        for ci_c0 in ci[pass_i]:
            ci0 = ci[pass_i][ci_c0]
            cis0 = cis[pass_i][ci_c0]
            
            if 'template_key' not in ci0:       continue
            tk0 = ci0['template_key']['subtomogram']
            if tk0 not in tkd:  continue        # this means the cluster is not selected in current pass, ignore this cluster

            c0 = tkd[tk0]

            # collect the alignment scores of all the subtomograms of ci0
            ci0s = set( str(_['subtomogram']) for _ in ci0['data_json']  )      # set of subtomograms of ci0
            al0 = [_ for _ in al if _['vol_key']['subtomogram'] in ci0s]        # collection of alignment scores of subtomograms of ci0

            ss = {}     # scores collected for each template
            for c1 in tk:
                tk1 = tk[c1]['subtomogram']

                if tk1 not in tkd:  continue        # this means the template is already deleted because its corresponding cluster is not specific

                ss[c1] = N.array( [_['align'][c1]['score'] for _ in al0]  )

            for c1 in ss:
                if c1 is c0:        continue
                ind_t = N.logical_and(N.isfinite(ss[c0]), N.isfinite(ss[c1]))
                if ind_t.sum() < test_sample_num_min:    continue
                t_, p_ = SCS.wilcoxon(ss[c1][ind_t], ss[c0][ind_t])
                st[c0][c1] = {'median0':N.median(ss[c0][ind_t]), 'median1':N.median(ss[c1][ind_t]), 't':t_, 'p':p_}

    return st

if __name__ == '__main__':
    file_stat_file = '/home/rcf-47/mxu/ln/frequent_structure/data/out/method/imputation-classification/multi/classification/30/0.005/0221/out/file_stat.json'
    pass_i = 1

    main()

