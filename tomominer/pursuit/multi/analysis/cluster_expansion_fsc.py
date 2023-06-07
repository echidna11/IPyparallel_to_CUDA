#!/usr/bin/env python

import os
import sys
import pickle
import json
import time
import uuid
from collections import defaultdict
from multiprocessing import Pool

import classify.util as CU


# Greedly expand selected clusters to maximize FSC sum
# parameters:       ci: cluster_info        dj: data_json
def cluster_expansion_fsc(self, ci, dj, out_dir, exp_op, avg_op, cluster_parallel_processing=False):

    current_pass_i = len(ci) - 1

    start_time = time.time()

    # collect the records that maximally match each template
    st = defaultdict(list)
    for r in dj:
        st[ str(r['template']['subtomogram']) ].append(r)

    
    expanded = False

    if not cluster_parallel_processing:
        for i in ci:
            for c in ci[i]:
                tk = str(ci[i][c]['template_key']['subtomogram'])
                if tk not in st:  continue

                pr = {}
                cluster_id_max = max(_ for _ in ci[current_pass_i])
                pr['pass_i'] = current_pass_i
                pr['cluster'] = cluster_id_max + 1

                pr = cluster_expansion_fsc__single_template(self=self, ci=ci[i][c], dj=st[tk], ci_re=pr, out_dir=out_dir, exp_op=exp_op, avg_op=avg_op)
                if pr is None:  continue

                ci[pr['pass_i']][pr['cluster']] = pr

                expanded = True
    else:

        pool = Pool()

        pool_results = []
        for i in ci:
            for c in ci[i]:
                tk = str(ci[i][c]['template_key']['subtomogram'])
                if tk not in st:  continue

                pr = {}
                cluster_id_max = max(_ for _ in ci[current_pass_i])
                pr['pass_i'] = current_pass_i
                pr['cluster'] = cluster_id_max + 1

                pool_results.append( pool.apply_async(func=cluster_expansion_fsc__single_template, kwds={'self':self, 'ci':ci[i][c], 'dj':st[tk], 'ci_re':pr, 'out_dir':out_dir, 'exp_op':exp_op, 'avg_op':avg_op}) )


        for r in pool_results:
            pr = r.get()
            if pr is None:  continue

            ci[pr['pass_i']][pr['cluster']] = pr
            expanded = True


    print 'cluster_expansion_fsc() took time ', time.time() - start_time

    return expanded




# for each cluster average, order all best matched subtomograms (denoted as S) according to alignment score
# include into cluster subtomograms in S one by one, and see at which stage the FSCc gets maximum
# then make a new average by adding a random uuid key
# parameters:       ci: cluster_info for a particular cluster        dj: data_json for a particular json
def cluster_expansion_fsc__single_template(self, ci, dj, ci_re=None, exp_op=None, avg_op=None, out_dir=None):


    if exp_op is None:      exp_op = {}


    cluster_size_org = len(ci['data_json'])

    djs = sorted( dj, key=lambda _ : (- float(_['score']) ) )


    if False:
        # we seprate all matched subtomograms into two sets according to given cluster, each set is ordered according to the alignment score

        subtomograms_t = set( str(_['subtomogram']) for _ in ci['data_json'] )

        djs0 = [_ for _ in djs if _['subtomogram'] in subtomograms_t];      assert( len(djs0) == len(subtomograms_t) )          # this is the set of subtomograms that already inside the cluster in ci
        djs1 = [_ for _ in djs if _['subtomogram'] not in subtomograms_t]       # this are subtomograms that are not in ci

        if len(djs1) == 0:      return None

        djsc = djs0 + djs1

        exp_op['minimum_cluster_size'] = cluster_size_org         # in the case, we only choose extended clusters

    else:
        
        # evaluate in sequence according only to alignment scores
        djsc = djs


    if 'extension_ratio' in exp_op:        djsc = djsc[: min( len(djsc), int( round(cluster_size_org) * float(exp_op['extension_ratio']) ) )]       # we limit the amount of extension in order to constrain computation cost



    import tomominer.statistics.ssnr as SS
    ss = SS.ssnr_sequential(self=self, data_json=djsc)

    # find the best top number of subtomograms
    best = {}
    for i, fsc in enumerate(ss['fsc']):
        if 'minimum_cluster_size' in exp_op:
            if i < (exp_op['minimum_cluster_size'] - 1):   continue

        fscs = sum(fsc)
        #print i, fscs, '\t',

        if 'fscs' not in best:
            best['fscs'] = fscs

        if fscs >= best['fscs']:
            best['fscs'] = fscs
            best['i'] = i
            best['fsc'] = fsc

    if 'version' not in ci:
        ci['version'] = 0


    print 'cluster_expansion_fsc__single_template() ', 'original cluster size ', cluster_size_org, 'fsc', sum(ci['fsc']), 
    print 're-aligned', 
    if len(ss['fsc']) >= cluster_size_org:      print 'fsc', sum(ss['fsc'][cluster_size_org - 1]), 
    print 'expand to size ', best['i'], 'fsc', best['fscs']



    ci_re['fsc'] = best['fsc']

    ci_re['subtomograms'] = set( str(_['subtomogram']) for _ in djsc[:best['i']]   )
    ci_re['data_json'] = djsc[:best['i']]

    
    ci_re['template_key'] = {}
    ci_re['template_key']['subtomogram'] = os.path.join( out_dir, 'clus_vol_avg_%03d.mrc'%(ci_re['cluster']) )
    ci_re['template_key']['mask'] = os.path.join( out_dir, 'clus_mask_avg_%03d.mrc'%(ci_re['cluster']) )


    # create and save cluster average
    CU.vol_avg__global(self=self, data_json=ci_re['data_json'], template_key=ci_re['template_key'], op=avg_op)

    return ci_re


   

if __name__ == '__main__':

    
    op_file = sys.argv[1]
    with open(op_file) as f:        op = json.load(f)

    # walk through every subdir, find and record last passes
    selected_folders = []
    for root, sub_folders, files in os.walk(os.getcwd()):

        max_t = {}
        max_t['pass_i'] = -1
        for sub_folder in sub_folders:
            sub_folder_path = os.path.join(root, sub_folder)

            data_json_file = os.path.join(sub_folder_path, 'data_config.json')
            if not os.path.exists(data_json_file):    continue
 
            cluster_info_file = os.path.join(sub_folder_path, 'cluster_info.pickle')
            if not os.path.exists(cluster_info_file):    continue
    
           
            pass_i = int(sub_folder.split('_')[1])
            if pass_i > max_t['pass_i']:
                max_t['pass_i'] = pass_i
                max_t['subfolder_path'] = sub_folder_path
                max_t['data_json_file'] = data_json_file
                max_t['cluster_info_file'] = cluster_info_file

        if max_t['pass_i'] >= 0:        selected_folders.append(max_t)

    
    def process_one_pass(self, fnames, op, cluster_parallel_processing=False):

        # load files from one pass and process
        cache = Cache(cache_dir=os.getenv('CACHE_DIR'))
        
        with open(fnames['cluster_info_file']) as f:       ci = pickle.load(f)
        with open(fnames['data_json_file']) as f:   dj = json.load(f)

        cluster_expansion_fsc(self=self, ci=ci, dj=dj, avg_op=op, cluster_parallel_processing=cluster_parallel_processing)

        #with open( os.path.join(fnames['root'], 'cluster_info_extended.pickle'), 'wb') as f:        pickle.dump(ci, f)


    for max_t in selected_folders:
        print max_t
        process_one_pass(fnames=max_t, op=op['cluster']['averaging'], cluster_parallel_processing=True)



