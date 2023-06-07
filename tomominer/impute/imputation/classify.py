
import time
import uuid
import os
import logging
from collections import defaultdict
from pprint import pprint
import json

from classify_config import config_options, parse_data
from queue_master import QueueMaster
import tomo

from tomominer.cluster import kmeans

def load_rotate_impute(data):
    

# The classification pipeline using a python coordinated Master/Worker
# paradigm.
#
# This code is for the master, which appends tasks to a running work queue, and
# recieves back results.

def classify_main(opt_file, data_file, tmp_dir, qhost, qport):
    """
    """

    start_time = time.time()
    opt  = config_options()
    opt.parse_config(opt_file)

    data = parse_data(data_file)


    # determine the shape of the data.
    v = tomo.parse_mrc(data.itervalues())       # get first entry subtomogram to load
    vol_shape = v.shape
    print "Finished parsing and loading: %2.6f sec" % (time.time() - start_time)

    runner = QueueMaster(qhost, qport)

    for pass_i in range(1, opt.pass_num):

        print "Beginning pass #%d" % (pass_i)
        pass_start_time = time.time()

        pass_dir  = tmp_dir + '/' + '/pass_%03d' % (pass_i)

        # TODO: do something smarter if directory already exists, or cannot be made.
        try:
            os.mkdir(pass_dir)
        except:
            pass


        # mxu: todo: kmeans clustering with interpolation, store averages to local disk
        # first load and impute all data
        # then perform kmeans clustering
        # then get averaged clusters

        start_time = time.time()
        kmeans_res = kmeans(data, opt.cluster_kmeans_k)

        print "k-means: %2.6f sec" % (time.time() - start_time)






        # align against template of largest cluster

        tasks = []
        template_keys = {}
        for c,local_vm in cluster_global.items():
            vol_avg_out_key  = pass_dir + '/' + 'template_%d_%d_avg.mrc' %(pass_i, c,)
            mask_avg_out_key = pass_dir + '/' + 'template_%d_%d_mask.mrc' %(pass_i, c,)
            template_keys[c] = (vol_avg_out_key, mask_avg_out_key)

            if opt.cluster_use_fft_avg:
                tasks.append(runner.task('template_avg_fft_global', vol_shape, local_vm, vol_avg_out_key, mask_avg_out_key, True, opt.template_mask_ratio_cutoff))
            else:
                tasks.append(runner.task('template_avg_global', vol_shape, local_vm, cluster_sizes[c], vol_avg_out_key, mask_avg_out_key, True, opt.template_mask_ratio_cutoff))

        wait = [_ for _ in runner.run(tasks)]
        print "Cluster templates %2.6f sec" % (time.time() - start_time)

        #--------------------------------
        # align all averages to average of largest cluster.
        largest_cluster = max( [ (len(clusters[c]), c) for c in clusters ] )[1]

        v1,m1 = template_keys[largest_cluster]

        # Then do the translation/rotation and save all the templates to disk.
        largest_vol_key  = pass_dir + '/template_%d_%d_vol_aligned.mrc'  % (pass_i, largest_cluster)
        largest_mask_key = pass_dir + '/template_%d_%d_mask_aligned.mrc' % (pass_i, largest_cluster)


        tomo.write_mrc(tomo.parse_mrc(v1), largest_vol_key)
        tomo.write_mrc(tomo.parse_mrc(m1), largest_mask_key)

        start_time = time.time()
        tasks = []

        aligned_templates = {}

        for c, (v,m) in template_keys.items():
            if c is not largest_cluster:
                vol_aligned_out_key  = pass_dir + '/template_%d_%d_vol_aligned.mrc'  % (pass_i, c)
                mask_aligned_out_key = pass_dir + '/template_%d_%d_mask_aligned.mrc' % (pass_i, c)
                tasks.append(runner.task('combined_search', largest_vol_key, largest_mask_key, v, m, opt.L, vol_aligned_out_key, mask_aligned_out_key))
                aligned_templates[c] = (vol_aligned_out_key, mask_aligned_out_key)

        wait = [_ for _ in runner.run(tasks)]
        print "Align cluster templates. %2.6f sec" % (time.time() - start_time)

        aligned_templates = aligned_templates.values()

        tasks = []
        for vmal in data:

            # send each volume alignment as its own work unit.
            tasks.append(runner.task("align_to_templates", data[vmal], aligned_templates, opt.L))

        start_time = time.time()
        align_to_template_record = {}     # mxu: for each subtomogram record the best aligned template (id) together with align score, angle and displacement
        for res in runner.run(tasks):
            assert res.result.subtomo_key in data
            assert data[res.result.subtomo_key]['subtomogram'] is res.result.subtomo_key

            data[res.result.subtomo_key]['ang'] = res.result.ang
            data[res.result.subtomo_key]['loc'] = res.result.loc
            data[res.result.subtomo_key]['score'] = res.result.score
            data[res.result.subtomo_key]['template'] = res.result.tem_key
            data[res.result.subtomo_key]['template_mask'] = res.result.tem_msk_key

        print "Align all volumes to cluster_templates. %2.6f sec" % (time.time() - start_time)


        print "Entire pass: %2.6f sec" % (time.time() - pass_start_time)
        print "---------------------------"



if __name__ == '__main__':

    import sys
    qhost = sys.argv[1]
    qport = 5011

    param = sys.argv[2]
    data  = sys.argv[3]
    tmp_dir = sys.argv[4]


    classify_main(param, data, tmp_dir, qhost, qport)



'''
todo: 
add an option: 1) use interpolation for classification or averaging, 2) use masked alignment 
add an option: if kmeans-k is 1, then do averaging only without clustering
write code for fast alignment without masking
'''
