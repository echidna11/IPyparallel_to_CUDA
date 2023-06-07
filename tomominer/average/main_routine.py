#!/usr/bin/env python

import os, time, random, json, copy
import cPickle as pickle

import numpy as N

import tomominer.pursuit.multi.util as CU

import tomominer.io.file as IV




"""

simple iterative averaging and alignment, adapted from tomominer.pursuit/multi/run.py

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


    old_fsc_sum = float('-inf')
    for pass_i in range(op['options']['pass_num']):

        pass_dir_next = os.path.join(out_dir, 'pass_%03d' % (pass_i+1))
        if os.path.exists(os.path.join(pass_dir_next, 'cluster_info.pickle')):     continue

        
        print '-------------------------------------------------------'
        print 'Beginning pass #%d' % (pass_i)
        pass_start_time = time.time()


        pass_dir  = os.path.join(out_dir, 'pass_%03d' % (pass_i))
        if not os.path.exists(pass_dir) :            os.mkdir(pass_dir)


        data_json_file =  os.path.join(pass_dir, 'data.json')
        cluster__kmeans_sfsc__file = os.path.join(pass_dir, 'cluster__kmeans_sfsc.pickle')
        cluster_averaging_file = os.path.join(pass_dir, 'cluster_averaging.pickle')

        if os.path.isfile(data_json_file):

            print 'loading ' + data_json_file
            with open(data_json_file) as f:     data_json = json.load(f)

            if op['fsc_stop_check']:
                with open(cluster__kmeans_sfsc__file, 'rb') as f:   cluster__kmeans_sfsc__file__re = pickle.load(f)
                cluster_ssnr_fsc = cluster__kmeans_sfsc__file__re['cluster_ssnr_fsc']
                old_fsc_sum = cluster_ssnr_fsc['fsc'][0].sum()

            with open(cluster_averaging_file, 'rb') as f:    ca_re = pickle.load(f)
            template_keys = ca_re['template_keys']


            print '....  and go to next pass'

            continue


        cluster_ssnr_fsc = None
        if op['fsc_stop_check']:
            # calculate the FSC of the global average
            if os.path.exists(cluster__kmeans_sfsc__file):
                print 'loading', cluster__kmeans_sfsc__file
                with open(cluster__kmeans_sfsc__file, 'rb') as f:   cluster__kmeans_sfsc__file__re = pickle.load(f)
                cluster_ssnr_fsc = cluster__kmeans_sfsc__file__re['cluster_ssnr_fsc']
            else:
                cluster_ssnr_fsc = CU.cluster_ssnr_fsc(self=self, clusters={0:data_json}, n_chunk=n_chunk)
                with open(cluster__kmeans_sfsc__file, 'wb') as f:     pickle.dump({'cluster_ssnr_fsc':cluster_ssnr_fsc}, f, protocol=-1)

            # stop if there is no increase in FSC
            fsc_sum = cluster_ssnr_fsc['fsc'][0].sum()
            print 'FSC sum', fsc_sum

            if fsc_sum > old_fsc_sum:
                old_fsc_sum = fsc_sum
            else:
                print 'no improvement, stopping'
                break



        #--------------------------------------------------------------------
        # calculate global average
        if os.path.exists(cluster_averaging_file):
            print 'loading ' + cluster_averaging_file
            with open(cluster_averaging_file, 'rb') as f:    ca_re = pickle.load(f)
        else:

            cluster_ssnr_fsc = None
            if ('smooth' in op['cluster']['averaging']) and (cluster_ssnr_fsc is None):             cluster_ssnr_fsc = CU.cluster_ssnr_fsc(self=self, clusters={0:data_json}, n_chunk=n_chunk)

            ca_op = copy.deepcopy(op['cluster']['averaging'])
            ca_op['n_chunk'] = n_chunk
            ca_op['out_dir'] = os.path.join(pass_dir, 'clus_avg')

            if 'smooth' in op['cluster']['averaging']:        ca_op['smooth']['fsc'] = {_:cluster_ssnr_fsc['fsc'][_] for _ in cluster_ssnr_fsc['fsc']}

            ca_re = CU.cluster_averaging(self=self, clusters={0:data_json}, op=ca_op)
            ca_re['cluster_ssnr_fsc'] = cluster_ssnr_fsc

            with open(cluster_averaging_file, 'wb') as f:    pickle.dump(ca_re, f, protocol=1)

        template_keys = ca_re['template_keys']
        #print template_keys
        if len(template_keys) == 0:             return      None



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
            for rec in data_json:   tasks.append(self.runner.task( priority=task_priority, module='tomominer.pursuit.multi.util', method='align_to_templates', kwargs={'rec':rec, 'tem_keys':ca_re['template_keys'], 'align_op':op['align'], 'multiprocessing':False} ))

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

            # numpy arrays are not pickle-able.  Convert to a list of numbers.
            rec['angle'] = [_ for _ in res['align'][0]['angle']]
            rec['loc'] = [_ for _ in res['align'][0]['loc']]
            rec['score'] = res['align'][0]['score']
        
            data_json_new.append(rec)



        data_json = data_json_new
        with open(data_json_file, 'w') as f:    json.dump(data_json, f, indent=2)


        print "Entire pass: %2.6f sec" % (time.time() - pass_start_time)



    return template_keys[0]     # return the best global average



