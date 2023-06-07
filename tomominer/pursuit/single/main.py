'''
given a set of subtomograms, with current alignment,
iteratively determine the best set using different methods (sequential, hierarchical clustering etc)
then re-average and re-align against the average
perform template guided segmentation if necessary

'''



import os, sys, json, copy, time, uuid, shutil
import cPickle as pickle

import numpy as N
from collections import defaultdict

import tomominer.pursuit.multi.util as PMU
import tomominer.pursuit.single.sequential.ssnr as PSSS


def stop_test_stat(op, fsc_stat, pass_i):
    print 'stop_test_stat()'

    if pass_i < 1:    return {'should_stop':False}
    
    best = None
    for i in fsc_stat:
        if i < 0:       continue
        if (best is None) or (best['fsc'] < fsc_stat[i]['fsc']):            best = fsc_stat[i] 

    assert  best is not None
    max_pass_i = best['pass_i']
    print 'max_pass_i', max_pass_i

    if (pass_i - max_pass_i) >= op['cluster']['stopping_test_pass_num']:
        return {'should_stop':True, 'max_pass_i':max_pass_i}
    else:
        return {'should_stop':False, 'max_pass_i':max_pass_i}


def pursuit(self, op, data_json):

    out_dir = os.path.abspath(op['out_dir'])
    if not os.path.isdir(out_dir):      os.makedirs(out_dir)

    print 'pursuit()', 'out_dir', out_dir

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
        # if no initial angles are given, set random ones
        for r in data_json:
            if 'loc' not in r: r['loc'] = [0.0, 0.0, 0.0]
            if 'angle' not in r: r['angle'] = [_ for _ in N.random.random(3) * (N.pi * 2)]
        del r
        with open(data_json_file, 'wb') as f:       json.dump(data_json, f, indent=2)


    n_chunk=max(op['options']['min_chunk_size'], int( N.ceil(len(data_json) / op['options']['worker_num']) ))

    file_stat_file = os.path.join(out_dir, 'file_stat.json')
    if os.path.isfile(file_stat_file):
        with open(file_stat_file) as f:        file_stat = json.load(f)
        file_stat['passes'] = {int(_):file_stat['passes'][_] for _ in file_stat['passes']}

    else:

        file_stat = {}
        file_stat['out_dir'] = out_dir
        file_stat['file_stat_file'] = file_stat_file
        file_stat['data_json_file'] = data_json_file     # this is the initial data_json_file
        file_stat['passes'] = {}
        file_stat['pass_i_current'] = -1
        with open(file_stat_file, 'w') as f:        json.dump(file_stat, f, indent=2)




    cluster_info = defaultdict(dict)       # use this to collect all cluster averages and their FSC information for cluster selection purpose. IMPORTANT: we do not store those information that are subject to change to here, for example 'is_specific'

    fsc_stat = {}
    fsc_stat[-1] = {'template_key':     op['initial_template']}       # the initial template key and fsc score


    for pass_i in range(op['options']['pass_num']):

        if pass_i < file_stat['pass_i_current']:   continue        # we only need data from the last completed pass
        
        print '-------------------------------------------------------'
        print 'Beginning pass #%d' % (pass_i)       ;       sys.stdout.flush()
        pass_start_time = time.time()

        pass_dir  = os.path.join(out_dir, 'pass', '%03d' % (pass_i))
        if not os.path.exists(pass_dir) :            os.makedirs(pass_dir)

        data_json_file =  os.path.join(pass_dir, 'data.json')
        fsc_stat_file = os.path.join(pass_dir, 'fsc_stat.json')

        if os.path.isfile(fsc_stat_file):           # this means that this pass is finished, load all data to prepare for the next pass

            # first load fsc_stat to determine whether need to stop
            print 'loading', fsc_stat_file
            with open(fsc_stat_file, 'r') as f:    fsc_stat = json.load(f)
            fsc_stat = fsc_stat_json_convert(fsc_stat)

            sts_re = stop_test_stat(op=op, fsc_stat=fsc_stat, pass_i=pass_i)

            if sts_re['should_stop']:
                print 'no more improvements seen, stop'
                break


            # load other files
            with open(data_json_file) as f:     data_json = json.load(f)
            print '....  and go to next pass'

            continue



        #assert      pass_i not in file_stat['passes']
        file_stat['passes'][pass_i] = {}
        file_stat['passes'][pass_i]['pass_i'] = pass_i
        file_stat['passes'][pass_i]['pass_dir'] = pass_dir
        file_stat['pass_i_current'] = pass_i


        #assert      pass_i not in fsc_stat
        fsc_stat[pass_i] = {}
        fsc_stat[pass_i]['pass_i'] = pass_i


        #--------------------------------------------------------------
        # decide the optimal set of subtomograms for averaging, and calculate SSNR


        sequential_expansion_file = os.path.join(pass_dir, 'sequential_expansion_file.pickle')
        file_stat['passes'][pass_i]['sequential_expansion_file'] = sequential_expansion_file

        if os.path.exists(sequential_expansion_file):
            print 'loading ', sequential_expansion_file
            with open(sequential_expansion_file, 'rb') as f:    tmp = pickle.load(f)
            sem = tmp['sem']
            del tmp
        else:
            sem_op = copy.deepcopy(op['cluster']['sequential_expansion'])
            sem_op['ssnr'] = copy.deepcopy(op['cluster']['ssnr'])
            if sem_op['ssnr']['segmentation'] == True:
                sem_op['ssnr']['segmentation'] = copy.deepcopy(op['segmentation'])
                sem_op['ssnr']['segmentation']['density_positive'] = op['density_positive']
            sem_op['averaging'] = copy.deepcopy(op['cluster']['averaging'])
            sem_op['min_expansion_size'] = op['cluster']['size_min']

            sem_op['n_chunk'] = n_chunk
            if False:
                sem_op['template'] = {'key':fsc_stat[pass_i-1]['template_key'], 'fsc':fsc_stat[pass_i-1]['fsc']}
            else:
                sem_op['template'] = {'key':fsc_stat[pass_i-1]['template_key']}

            
            if sem_op['mode'] == 'local_maxima':
                sem = PSSS.ssnr_sequential_expansion__local_maxima(self=self, dj=data_json, op=sem_op)
            elif sem_op['mode'] == 'global_maximum':
                sem = PSSS.ssnr_sequential_expansion__global_maxia(self=self, dj=data_json, op=sem_op)

            with open(sequential_expansion_file, 'wb') as f:    pickle.dump({'sem_op':sem_op, 'sem':sem}, f, protocol=1)
        
        fsc_stat[pass_i]['fsc'] = sem['fsc']

        print 'resulting cluster size', len(sem['data_json']), 'fsc', sem['fsc']

        #--------------------------------------------------------------------
        # Using the set we just generated, we build averages.

        cluster_averaging_file = os.path.join(pass_dir, 'cluster_averaging.pickle')
        file_stat['passes'][pass_i]['cluster_averaging_file'] = cluster_averaging_file

        if os.path.exists(cluster_averaging_file):
            print 'loading ', cluster_averaging_file
            with open(cluster_averaging_file, 'rb') as f:    ca_re = pickle.load(f)
        else:
            print 'build cluster averages'
            cu_ca_op = copy.deepcopy(op['cluster']['averaging'])
            cu_ca_op['out_dir'] = os.path.join(pass_dir, 'clus_avg')
            cu_ca_op['pass_i'] = pass_i         # this is used to be recorded into the template_key
            cu_ca_op['n_chunk'] = n_chunk
            ca_re = PMU.cluster_averaging(self=self, clusters={0:sem['data_json']}, op=cu_ca_op)

            # segment those templates aligned to the common frame
            if segmentation_op is not None:
                print 'segmenting the selected and aligned templates'
                template_segmentation_op = copy.deepcopy(segmentation_op)
                if ('segmentation' in op['template']) and ('normalize_and_take_abs' in op['template']['segmentation']) and (op['template']['segmentation']['normalize_and_take_abs']):      template_segmentation_op['normalize_and_take_abs'] = True
                PMU.template_segmentation(self=self, tk=ca_re['template_keys'], op=template_segmentation_op)

            with open(cluster_averaging_file, 'wb') as f:    pickle.dump(ca_re, f, protocol=1)


        cluster_info[pass_i][0] = {}
        cluster_info[pass_i][0]['template_key'] = ca_re['template_keys'][0]


        fsc_stat[pass_i]['template_key'] = ca_re['template_keys'][0]
        with open(fsc_stat_file, 'w') as f:    json.dump(fsc_stat, f, indent=2)


        sys.stdout.flush()




        #----------------------------------------------------------
        # align against selected templates
        align_template_file = os.path.join(pass_dir, 'align_template.pickle')
        file_stat['passes'][pass_i]['align_template_file'] = align_template_file

        if os.path.isfile(align_template_file):
            print 'loading ' + align_template_file
            with open(align_template_file, 'rb') as f:    at_ress = pickle.load(f)

        else:

            align_template__tmp_dir__file = os.path.join(pass_dir, 'align_template__tmp_dir.json')

            if os.path.isfile(align_template__tmp_dir__file):
                with open(align_template__tmp_dir__file) as f:       align_template__tmp_dir = json.load(f)
            else:
                align_template__tmp_dir = os.path.join(self.cache.tmp_dir, 'align-template-'+str(uuid.uuid4()))
                with open(align_template__tmp_dir__file, 'w') as f:       json.dump(align_template__tmp_dir, f, indent=2)

            start_time = time.time()

            if not os.path.isdir(align_template__tmp_dir):  os.makedirs(align_template__tmp_dir)

            at_ress = PMU.align_to_templates__batch(self=self, op=op, data_json=data_json, segmentation_tg_op=segmentation_tg_op, tmp_dir=align_template__tmp_dir, tem_keys={0:cluster_info[pass_i][0]['template_key']})

            with open(align_template_file, 'wb') as f:    pickle.dump(at_ress, f, protocol=1)
            shutil.rmtree(align_template__tmp_dir)          #  clean up temporary alignment info

            print "Align all volumes to cluster_templates. %2.6f sec" % (time.time() - start_time)

        at_ress = [_.result for _ in at_ress]


        sys.stdout.flush()

        cluster_info_file = os.path.join(pass_dir, 'cluster_info.pickle')
        file_stat['passes'][pass_i]['cluster_info_file'] = cluster_info_file
        with open(cluster_info_file, 'wb') as f:    pickle.dump(cluster_info[pass_i], f, protocol=-1)           # we only dump cluster info of CURRENT PASS, to save storage!!

        # record template alignment results
        data_json_new = []
        for res in at_ress:

            r = {}
            r['subtomogram'] = res['vol_key']['subtomogram']
            r['mask'] = res['vol_key']['mask']

            # numpy arrays are not pickle-able.  Convert to a list of numbers.
            r['angle'] = [_ for _ in res['align'][0]['angle']]
            r['loc'] = [_ for _ in res['align'][0]['loc']]
            r['score'] = res['align'][0]['score']
            
            data_json_new.append(r)
            del r


        data_json = data_json_new
        with open(data_json_file, 'w') as f:    json.dump(data_json, f, indent=2)
        file_stat['passes'][pass_i]['data_json_file'] = data_json_file

        sts_re = stop_test_stat(op=op, fsc_stat=fsc_stat, pass_i=pass_i)
        if 'pass_i_max' in file_stat:           file_stat['pass_i_max'] = sts_re['max_pass_i']

        with open(file_stat_file, 'w') as f:        json.dump(file_stat, f, indent=2)


        if sts_re['should_stop']:
            print 'no more improvements seen, stop'
            break


    return file_stat

