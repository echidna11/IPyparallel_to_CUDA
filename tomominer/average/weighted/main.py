#!/usr/bin/env python

import os, time, random, json, copy
import cPickle as pickle

import numpy as N


import tomominer.io.file as IV

import tomominer.average.weighted.util as AWU



"""

simple iterative averaging and alignment, adapted from tomominer.pursuit/multi/run.py

we use weighted averaging instead of directly perform averaging, see in such case if we can reduce missing wedge problem...

in such case, we can expect that the SNR is largely reduced, because for each subtomogram we contribute multiple copies into the average and the noise pixels are rotated and averaged out... Therefore we probably do not need to perform low pass filtering to the average



~/ln/tomominer/tomominer/average/weighted/main.py

"""





def do_average(self, op, data_json, out_dir):

    out_dir = os.path.abspath(out_dir)
    if not os.path.isdir(out_dir):      os.makedirs(out_dir)

    if 'test' in op:
        if ('sample_num' in op['test']) and (op['test']['sample_num'] > 0) and (len(data_json) > op['test']['sample_num']):
            print 'testing the procedure using a subsample of %d subtomograms'%(op['test']['sample_num'])
            data_json = random.sample(data_json, op['test']['sample_num'])


    data_json_file = os.path.join(out_dir, 'data.json')
    if os.path.exists(data_json_file):
        print 'loading ' + data_json_file
        with open(data_json_file) as f:     data_json = json.load(f)

    else:
        # if no initial angles are given, set random ones
        for rec in data_json:
            if 'loc' not in rec: rec['loc'] = [0.0, 0.0, 0.0]
            if 'angle' not in rec: rec['angle'] = [_ for _ in N.random.random(3) * (N.pi * 2)]
        with open(data_json_file, 'wb') as f:       json.dump(data_json, f, indent=2)

    # get size of subtomogram
    v = IV.get_mrc(   data_json[0]['subtomogram']   )
    size = N.array(v.shape)
    del v


    n_chunk=max(op['options']['min_chunk_size'], int( N.ceil(len(data_json) / op['options']['worker_num']) ))
    n_chunk_quick_job=max(op['options']['min_chunk_size'], int( N.ceil(len(data_json) / op['options']['worker_num_quick_job']) ))       # use this to reduce those quick jobs that resulting large IO consumption


    for pass_i in range(op['options']['pass_num']):

        pass_dir_next = os.path.join(out_dir, 'pass_%03d' % (pass_i+1))
        if os.path.exists(os.path.join(pass_dir_next, 'cluster_info.pickle')):     continue

        
        print '-------------------------------------------------------'
        print 'Beginning pass #%d' % (pass_i)
        pass_start_time = time.time()


        pass_dir  = os.path.join(out_dir, 'pass_%03d' % (pass_i))
        if not os.path.exists(pass_dir) :            os.mkdir(pass_dir)


        data_json_file =  os.path.join(pass_dir, 'data.json')
        cluster_averaging_file = os.path.join(pass_dir, 'cluster_averaging.pickle')

        if os.path.isfile(data_json_file):

            print 'loading ' + data_json_file
            with open(data_json_file) as f:     data_json = json.load(f)

            with open(cluster_averaging_file, 'rb') as f:    avg_re = pickle.load(f)

            print '....  and go to next pass'

            continue



        #--------------------------------------------------------------------
        # calculate global average
        if os.path.exists(cluster_averaging_file):
            print 'loading ' + cluster_averaging_file
            with open(cluster_averaging_file, 'rb') as f:    avg_re = pickle.load(f)

        else:

            avg_op = copy.deepcopy(op['cluster']['averaging'])
            avg_op['n_chunk'] = n_chunk
            avg_op['out_dir'] = os.path.join(pass_dir, 'clus_avg')

            avg_re = AWU.cluster_averaging(self=self, clusters={0:data_json}, op=avg_op)

            with open(cluster_averaging_file, 'wb') as f:    pickle.dump(avg_re, f, protocol=1)



        #----------------------------------------------------------
        # align against selected templates
        align_template_file = os.path.join(pass_dir, 'align_template.pickle')
        if os.path.exists(align_template_file):
            print 'loading ' + align_template_file
            with open(align_template_file, 'rb') as f:    at_ress = pickle.load(f)
        else:
            print 'alignment', op['align']
            task_priority = 2000 + N.random.randint(100)        # add a random number so that this batch of tasks stay together
            tasks = []
            for rec in data_json:   tasks.append(self.runner.task( priority=task_priority, module='tomominer.average.weighted.util', method='align_to_templates', kwargs={'rec':rec, 'tem_keys':avg_re['template_keys'], 'align_op':op['align'], 'multiprocessing':False} ))

            start_time = time.time()
            at_ress = [_ for _ in self.runner.run__except(tasks)]

            with open(align_template_file, 'wb') as f:    pickle.dump(at_ress, f, protocol=1)

            print "Align all volumes to cluster_templates. %2.6f sec" % (time.time() - start_time)

        at_ress = [_.result for _ in at_ress]



        # record template alignment results
        data_json_new = []
        for res in at_ress:

            rec = {}
            rec['subtomogram'] = res['vol_key']['subtomogram']
            rec['mask'] = res['vol_key']['mask']

            # here we record all alignments
            rec['align'] = []
            for a in res['align'][0]['align']:
                at = {}
                at['score'] = a['score']
                at['loc'] = a['loc'].tolist()
                at['angle'] = a['angle'].tolist()

                rec['align'].append(at)
                del at
       
            data_json_new.append(rec)



        data_json = data_json_new
        with open(data_json_file, 'w') as f:    json.dump(data_json, f, indent=2)


        print "Entire pass: %2.6f sec" % (time.time() - pass_start_time)



    return avg_re['template_keys'][0]     # return the best global average



