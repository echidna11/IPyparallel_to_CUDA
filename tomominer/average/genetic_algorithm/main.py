#!/usr/bin/env python

import os, time, random, json, copy
import cPickle as pickle
import numpy as N
import tomominer.io.file as IV
import tomominer.pursuit.multi.util as PMU
import tomominer.pursuit.multi.pattern_generate.genetic_algorithm_ssnr_fsc as PMPG

"""
simple iterative averaging and alignment, adapted from tomominer.pursuit/multi/run.py
we use genetic algorithm to find a subset of subtomograms that achieve maximum quality / resolution, then align all subtomograms against the average of this subset
~/ln/tomominer/tomominer/average/genetic_algorithm/main.py
"""

def do_average(self, op, data_json):

    out_dir = os.path.abspath(op['out_dir'])
    if not os.path.isdir(out_dir):      os.makedirs(out_dir)

    if 'test' in op:
        if ('sample_num' in op['test']) and (op['test']['sample_num'] > 0) and (len(data_json) > op['test']['sample_num']):
            print 'testing the procedure using a subsample of %d subtomograms'%(op['test']['sample_num'])
            data_json = random.sample(data_json, op['test']['sample_num'])
    else:
        print 'total subtomogram number', len(data_json)



    if 'segmentation' in op:
        segmentation_op = copy.deepcopy(op['segmentation'])
        segmentation_op['density_positive'] = op['density_positive']

        segmentation_tg_op = copy.deepcopy(segmentation_op)     # options for segmenting a subtomogram, guided by template
        if ('guided_segmentation' in op['template']) and ('gaussian_smooth_sigma' in op['template']['guided_segmentation']):            segmentation_tg_op['gaussian_smooth_sigma'] =  op['template']['guided_segmentation']['gaussian_smooth_sigma']

    else:
        segmentation_op = None
        segmentation_tg_op = None



    data_json_file = os.path.join(out_dir, 'data.json')
    if os.path.exists(data_json_file):
        print 'loading ' + data_json_file
        with open(data_json_file) as f:     data_json = json.load(f)

    else:
        if op['random_initial_alignment']:
            print 'force use random initial orientations'
        else:
            print 'use initial transformations if avaliable'

        # if no initial angles are given, set random ones
        random_assign_count = 0
        for rec in data_json:
            randomly_assigned = False

            if ('loc' not in rec) or op['random_initial_alignment']: 
                rec['loc'] = [0.0, 0.0, 0.0]
                randomly_assigned = True

            if ('angle' not in rec) or op['random_initial_alignment']: 
                rec['angle'] = [_ for _ in N.random.random(3) * (N.pi * 2)]
                randomly_assigned = True

            if randomly_assigned:       random_assign_count += 1

        with open(data_json_file, 'wb') as f:       json.dump(data_json, f, indent=2)

        print 'randomly assigned', random_assign_count, 'orientations', 'out of', len(data_json), 'subtomograms'

    # get size of subtomogram
    v = IV.get_mrc(   data_json[0]['subtomogram']   )
    size = N.array(v.shape)
    del v

    file_stat_file = os.path.join(out_dir, 'file_stat.json')
    if os.path.isfile(file_stat_file):
        with open(file_stat_file) as f:        file_stat = json.load(f)
        file_stat['passes'] = {int(_):file_stat['passes'][_] for _ in file_stat['passes']}

    else:

        file_stat = {}
        file_stat['out_dir'] = out_dir
        file_stat['file_stat_file'] = file_stat_file
        file_stat['data_json_file'] = data_json_file     # this is the initial data_json_file
        file_stat['pass_i_current'] = 0
        file_stat['passes'] = {}
        with open(file_stat_file, 'w') as f:        json.dump(file_stat, f, indent=2)




    pmpg_old = None
    fsc_stat = {}

    for pass_i in range(op['options']['pass_num']):

        if pass_i < (file_stat['pass_i_current'] - 1):   continue        # we only need data from the last completed pass
 

        print '-------------------------------------------------------'
        print 'Beginning pass #%d' % (pass_i)
        pass_start_time = time.time()


        pass_dir  = os.path.join(out_dir, 'p', '%03d' % (pass_i))
        if not os.path.exists(pass_dir) :            os.makedirs(pass_dir)


        data_json_file =  os.path.join(pass_dir, 'data.json')
        pmpg_file = os.path.join(pass_dir, 'pmpg.pickle')
        fsc_stat_file = os.path.join(pass_dir, 'fsc_stat.json')
        cluster_averaging_file = os.path.join(pass_dir, 'cluster_averaging.pickle')

        if os.path.isfile(data_json_file):

            print 'loading', data_json_file
            with open(data_json_file) as f:     data_json = json.load(f)

            print 'loading', pmpg_file
            with open(pmpg_file, 'rb') as f:        pmpg_old = pickle.load(f)

            print 'loading', fsc_stat_file
            with open(fsc_stat_file) as f:     fsc_stat = json.load(f)
            fsc_stat = {int(_):fsc_stat[_] for _ in fsc_stat}

            print 'loading', cluster_averaging_file
            with open(cluster_averaging_file, 'rb') as f:    avg_re = pickle.load(f)


            print '....  and go to next pass'

            continue


        n_worker = max(self.runner.work_queue.get_worker_number(), 1)
        n_chunk = max(op['options']['min_chunk_size'], int( N.ceil(len(data_json) / n_worker) ))

        file_stat['pass_i_current'] = pass_i
        if pass_i not in file_stat['passes']:
            file_stat['passes'][pass_i] = {}
        with open(file_stat_file, 'w') as f:        json.dump(file_stat, f, indent=2)


        #--------------------------------------------------------------------
        # select best subset
        if not os.path.isfile(pmpg_file):
            pmpg_dp_op = {}
            if segmentation_tg_op is not None:  pmpg_dp_op['segmentation_tg'] = segmentation_tg_op


            pmpg_ga_op = copy.deepcopy(op['genetic_algorithm'])
            pmpg_ga_op['sum_min'] = op['min_sample_num']
            pmpg_ga_op['evaluate']['ssnr']['mask_sum_threshold'] = op['min_sample_num']

            if 'parallel' in pmpg_ga_op:        pmpg_ga_op['parallel']['n_chunk'] = n_chunk

            pmpg = {}
            pmpg['dp'] = PMPG.data_prepare(dj=data_json, op=pmpg_dp_op)
            pmpg['best'] = PMPG.ga(self=self, stat=pmpg['dp'], initial_population=(pmpg_old['best']['p'] if pmpg_old is not None else None), op=pmpg_ga_op)
            pmpg['dj'] = [data_json[_] for _ in range(len(data_json)) if (pmpg['best']['p'][0,_]==1)]

            del pmpg['dp']          # save some storage
            with open(pmpg_file, 'wb') as f:        pickle.dump(pmpg, f, protocol=-1)

        else:

            with open(pmpg_file, 'rb') as f:        pmpg = pickle.load(f)


        pmpg_old = pmpg

        if 'pmpg_file' not in file_stat['passes'][pass_i]:
            file_stat['passes'][pass_i]['pmpg_file'] = pmpg_file
            with open(file_stat_file, 'w') as f:        json.dump(file_stat, f, indent=2)


        if os.path.exists(fsc_stat_file):
            print 'loading', fsc_stat_file
            with open(fsc_stat_file) as f:     fsc_stat = json.load(f)
            fsc_stat = {int(_):fsc_stat[_] for _ in fsc_stat}
        else:
            fsc_stat[pass_i] = {}
            fsc_stat[pass_i]['fsc_sum'] = N.sum(    pmpg['best']['e'][0]['fsc']     )
            fsc_stat[pass_i]['cluster_size'] = pmpg['best']['p'][0,:].sum()

            with open(fsc_stat_file, 'w') as f:         json.dump(fsc_stat, f, indent=2)

        if 'fsc_stat_file' not in file_stat['passes'][pass_i]:
            file_stat['passes'][pass_i]['fsc_stat_file'] = fsc_stat_file
            with open(file_stat_file, 'w') as f:        json.dump(file_stat, f, indent=2)

        stopping_test_re = stopping_test(fs=fsc_stat, op=op['stop'])
        file_stat['best_pass_i'] = stopping_test_re['best_pass_i']
        with open(file_stat_file, 'w') as f:        json.dump(file_stat, f, indent=2)

        print 'best_pass_i', file_stat['best_pass_i'], 'FSC score', fsc_stat[file_stat['best_pass_i']]['fsc_sum'], 'cluster_size', fsc_stat[file_stat['best_pass_i']]['cluster_size']

        if stopping_test_re['should_stop']:
            print 'stopping'
            return file_stat

        #--------------------------------------------------------------------
        # calculate global average
        if os.path.exists(cluster_averaging_file):
            print 'loading ' + cluster_averaging_file
            with open(cluster_averaging_file, 'rb') as f:    avg_re = pickle.load(f)

        else:

            avg_op = copy.deepcopy(op['averaging'])
            avg_op['mask_count_threshold'] = op['min_sample_num']
            avg_op['n_chunk'] = n_chunk
            avg_op['out_dir'] = os.path.join(pass_dir, 'avg')
            if 'smooth' in avg_op:                avg_op['smooth']['fsc'] = {0:pmpg['best']['e'][0]['fsc']}

            avg_re = PMU.cluster_averaging(self=self, clusters={0:pmpg['dj']}, op=avg_op)

            with open(cluster_averaging_file, 'wb') as f:    pickle.dump(avg_re, f, protocol=1)

        if 'cluster_averaging_file' not in file_stat['passes'][pass_i]:
            file_stat['passes'][pass_i]['cluster_averaging_file'] = cluster_averaging_file
            with open(file_stat_file, 'w') as f:        json.dump(file_stat, f, indent=2)

        # segment the average
        if segmentation_op is not None:
            print 'segmenting the selected and aligned templates'
            template_segmentation_op = copy.deepcopy(segmentation_op)
            if ('segmentation' in op['template']) and ('normalize_and_take_abs' in op['template']['segmentation']) and (op['template']['segmentation']['normalize_and_take_abs']):      template_segmentation_op['normalize_and_take_abs'] = True
            PMU.template_segmentation(self=self, tk=avg_re['template_keys'], op=template_segmentation_op)


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
            for rec in data_json:   tasks.append(self.runner.task( priority=task_priority, module='tomominer.pursuit.multi.util', method='align_to_templates', kwargs={'rec':rec, 'tem_keys':avg_re['template_keys'], 'align_op':op['align'], 'segmentation_tg_op':(segmentation_tg_op if op['template']['match']['use_segmentation_mask'] else None), 'multiprocessing':False} ))

            start_time = time.time()
            at_ress = [_ for _ in self.runner.run__except(tasks)]

            with open(align_template_file, 'wb') as f:    pickle.dump(at_ress, f, protocol=1)

            print "Align all volumes to cluster_templates. %2.6f sec" % (time.time() - start_time)

        if 'align_template_file' not in file_stat['passes'][pass_i]:
            file_stat['passes'][pass_i]['align_template_file'] = align_template_file
            with open(file_stat_file, 'w') as f:        json.dump(file_stat, f, indent=2)


        at_ress = [_.result for _ in at_ress]



        # record template alignment results
        data_json_new = []
        for res in at_ress:

            rec = {}
            rec['subtomogram'] = res['vol_key']['subtomogram']
            rec['mask'] = res['vol_key']['mask']
            rec['score'] = res['align'][0]['score']
            rec['loc'] = res['align'][0]['loc'].tolist()
            rec['angle'] = res['align'][0]['angle'].tolist()

            data_json_new.append(rec)



        data_json = data_json_new
        with open(data_json_file, 'w') as f:    json.dump(data_json, f, indent=2)

        print 'failure alignment number', len([_ for _ in data_json if not N.isfinite(_['score'])])


        if 'data_json_file' not in file_stat['passes'][pass_i]: 
            file_stat['passes'][pass_i]['data_json_file'] = data_json_file
            with open(file_stat_file, 'w') as f:        json.dump(file_stat, f, indent=2)

        print "Entire pass: %2.6f sec" % (time.time() - pass_start_time)



    return file_stat


def stopping_test(fs, op):

    pass_i = N.max( [_ for _ in fs] )

    best_pass_i = 0
    for i in fs:
        if fs[i]['fsc_sum'] <= fs[best_pass_i]['fsc_sum']:   continue
        best_pass_i = i

    if (pass_i - best_pass_i) > op['max_none_increasing_iteration_num']:
        return {'should_stop':True, 'best_pass_i':best_pass_i}
    else:
        return {'should_stop':False, 'best_pass_i':best_pass_i}



