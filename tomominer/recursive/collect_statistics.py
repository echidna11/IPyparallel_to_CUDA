

# recursively walk through clustering result and collect all information for statistics, 
# store the result into a matlab file for analysis under matlab

import os
import sys

import json
import scipy.io as sio

import classify_config

from hierarchy_coding_file_name import extract_hierarchy_code
from subtomogram_cluster_sel import tomominer.cluster_sel
from work_queue.queue_master import QueueMaster

import vol

import tomo

if __name__ == '__main__':
    
    qhost = sys.argv[1]
    qport = 5011


    runner = QueueMaster(qhost, qport)



    op_file = sys.argv[2]       # option file for recursive classification, only used for extracting the pass number and cluster number

    src_dir = sys.argv[3]
    tmp_dir = sys.argv[4]
    res_dir = sys.argv[5]

    op = classify_config.config_options()
    op.parse_config(op_file)

    vol_shape = None

    for root, subFolders, files in os.walk(src_dir):

        pass_s = root.split('/')
        pass_s = pass_s[len(pass_s)-1]
        if not pass_s.startswith('pass_'):  continue
        pass_s = pass_s.split('_')    

        try:
            pass_i = int(pass_s[1])
        except ValueError:
            continue    
        
        if pass_i != (op.pass_num - 1)  :   continue        # we only process the last pass
            

        (hc_s, hc_i) = extract_hierarchy_code(root)
        if len(hc_i) == 0:        continue
        print 'hierarchy code ' + hc_s,     

        result_file = os.path.join(res_dir, hc_s + '.mat')
        if os.path.isfile(result_file):
            print 'ignoring ' + result_file
            continue

        data_json_file = None
        for f in files:
            if not f.endswith('.json'):     continue
            if not f.startswith('data_config_'): continue
            data_json_file = f
            break
        
        if data_json_file is None:  continue                

        # load json file that contains subtomogram alignment and clustering info
        with open(os.path.join(root, data_json_file)) as f:
            data_json = json.load(f)

        clus_info_s = []
        for clus_i in range(op.cluster_kmeans_k):
            print '  ' + repr(clus_i),

            clus_info = {}
            
            clus_info['hierarchy'] = hc_i
            clus_info['hierarchy_s'] = hc_s
            clus_info['cluster'] = clus_i

            data_json_t = cluster_sel(data_json, set([clus_i]), include=True)
             
            # calculate cluster size
            clus_info['size'] = len(data_json_t)

            # load and store cluster averages (and corresponding masks)
            avg_file = os.path.join(root, 'template_%d_%d_avg.mrc'%(pass_i, clus_i))
            if os.path.isfile(avg_file):
                clus_info['avg'] = tomo.parse_mrc(avg_file)


            avg_mask_file = os.path.join(root, 'template_%d_%d_mask.mrc'%(pass_i, clus_i))
            if os.path.isfile(avg_mask_file):
                clus_info['avg_maks'] = tomo.parse_mrc(avg_mask_file)

            data = classify_config.parse_data(data_json)

            if vol_shape is None:
                # determine the shape of the data.
                v = tomo.parse_mrc(data[0][0])
                vol_shape = v.shape

            # seperate a cluster into two halves, calculate average for each half
            ah = vol.average_halves({'data':data, 'runner':runner, 'tmp_dir':tmp_dir, 'vol_shape':vol_shape})

            ## do some format conversion, so that the data can be saved to matlab format
            data_h = ah['data_halves']
            data_h_m = []
            for data_h_t in data_h:
                data_h_m_t = []
                for v, m, a, l in data_h_t:
                    data_h_m_t.append({'vol':v, 'mask':m, 'ang':a, 'loc':l})
                data_h_m.append(data_h_m_t)    

            clus_info['average_halves'] = {'data_halves':data_h_m, 'averages':ah['averages']}

            clus_info_s.append(clus_info)

        print

        # store all collected info into a matlab mat file for further analysis, we keep saving for every iteration, because the collection process is slow
        sio.savemat(result_file, {'clus_info':clus_info_s})


