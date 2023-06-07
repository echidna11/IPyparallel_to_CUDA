

"""

pursuit, classification and averaging with pure missing wedge inerpolation,   no missing wedge correction

in this version, we use K-means clustering and sequential expansion etc to generate subtomogram sets.
In such case, the cluster SSNR are calculated simultaneously. 
~/ln/tomominer/tomominer/pursuit/multi/main.py
"""

import os, sys, shutil, time, copy, random, json, uuid
from collections import defaultdict
import cPickle as pickle
import numpy as N
import tomominer.pursuit.multi.util as CU
import tomominer.io.file as IV

# convert the keys from string to int
def fsc_stat_json_convert(fsc_stat):
    fsc_stat_t = {}
    for pass_i, fsc_stat_ in fsc_stat.iteritems():
        fsc_stat_t[int(pass_i)] = {}

        for c, fsc_stat__ in fsc_stat_.iteritems():
            fsc_stat_t[int(pass_i)][int(c)] = fsc_stat__

    return fsc_stat_t


def stop_test_stat(op, fsc_stat, pass_i, cluster_modes):
    #print 'stop_test_stat()'

    if op['options']['stopping_test']['criterion'] == 0:
        #print '# stopping test according to selected_templates_max_pass_i (of only non-redundant clusters)'
        #selected_templates_max_pass_i = N.max( [ cas_re['tk_info'][cas_re['selected_templates'][_]['subtomogram']]['pass_i'] for _ in cas_re['selected_templates'] ] )
        #selected_templates_max_pass_i = N.max(  [N.max( [__fsc_stat['pass_i'] for __, __fsc_stat in _fsc_stat.iteritems() if (__fsc_stat['is_specific'] is None)] ) for _, _fsc_stat in fsc_stat.iteritems()]  )      # we only consider specific clusters, in order to avoid the situation that new clusters populated then marked as non-specific then populated out again!!!           

        selected_templates_max_pass_i = 0
        for _, _fsc_stat in fsc_stat.iteritems():
            for  __, __fsc_stat in _fsc_stat.iteritems():
                if __fsc_stat['is_specific'] is not None:   continue        # we only consider specific clusters, in order to avoid the situation that new clusters populated then marked as non-specific then populated out again!!!
                if __fsc_stat['pass_i'] > selected_templates_max_pass_i:
                    selected_templates_max_pass_i = __fsc_stat['pass_i']
                    
        selected_templates_min_pass_i = N.min(  [_fsc_stat['pass_i'] for _, _fsc_stat in fsc_stat[pass_i].iteritems() if (_fsc_stat['is_specific'] is None)] )      # we only consider specific clusters selected in current pass 

        del _, _fsc_stat, __, __fsc_stat
        #print 'selected_templates_max_pass_i', selected_templates_max_pass_i
        print 'max pass of selected non-redundant average', selected_templates_max_pass_i

    else:               raise Exception('options_stopping_test_criterion')

    should_stop = False
    # This checks whether the last run of this program have already done 'stopping_test_pass_num'=5 iterations of kmeans without any improvement
    if (pass_i - selected_templates_max_pass_i) >= op['cluster']['stopping_test_pass_num']:        should_stop = True

    cluster_modes = copy.deepcopy(cluster_modes)
    print cluster_modes
    # This checks that if stopping_test_pass_num=5 iterations of kmeans are already done without improvement then whether or not adaptive_kmeans and sequential expansion have ben done or not.
    if should_stop:

        '''
        in order to prevent the premature of the clusters as much as possible, we first run kmeans mode as much as possible until no improvements under kmeans mode, then start other modes for refinements.
        '''

        if ('adaptive_k_ratio' in op['cluster']['kmeans']) and ('kmeans-adaptave' not in set(cluster_modes)):
            print 'include kmeans-adaptave mode for cluster population at next iteration'
            cluster_modes.append('kmeans-adaptave')
            should_stop = False

        if ('sequential_expansion' in op['cluster']) and ('sequential' not in set(cluster_modes)):
            print 'include sequential clustering mode for cluster population at next iteration'
            cluster_modes.append('sequential')
            should_stop = False

        if not should_stop:     print 'Starting next stage'

    return {'selected_templates_max_pass_i':selected_templates_max_pass_i, 'selected_templates_min_pass_i':selected_templates_min_pass_i, 'cluster_modes':cluster_modes, 'should_stop':should_stop}



def pursuit(self, op, data_json):
    # Create output directory if it doesn't exist
    out_dir = os.path.abspath(op['out_dir'])
    print 'Running pursuit()'
    print 'out_dir:', out_dir
    if not os.path.isdir(out_dir):      os.makedirs(out_dir)

    # For testing purpose, randomly sample the dataset of subtomograms and proceed. Sample number is mentiones in parameter file pursuit-op.json under tag "test"
    if 'test' in op:
        #if ('sample_num' in op['test']) and (op['test']['sample_num'] > 0) and (len(data_json) > op['test']['sample_num']):
        if ('sample_num' in op['test']) and (op['test']['sample_num'] > 0):
            print 'Testing the procedure using a subsample of %d subtomograms'%(op['test']['sample_num'])
            data_json = random.sample(data_json, op['test']['sample_num'])

    # Segmentation
    if 'segmentation' in op:
        if 'phi_propotion_cutoff' in op['segmentation']:
            raise Exception('The phi_propotion_cutoff option has moved to template.guided_segmentation for better organization')

        # Copy segmentation parameters
        segmentation_op = copy.deepcopy(op['segmentation'])
        segmentation_op['density_positive'] = op['density_positive']

        # Make another copy for template guided segmentation and copy parameters from tag "template">>"guided_segmentation" as well
        segmentation_tg_op = copy.deepcopy(segmentation_op)     # options for segmenting a subtomogram, guided by template
        if 'guided_segmentation' in op['template']:
            if 'gaussian_smooth_sigma' in op['template']['guided_segmentation']:                segmentation_tg_op['gaussian_smooth_sigma'] =  op['template']['guided_segmentation']['gaussian_smooth_sigma']
            if 'phi_propotion_cutoff' in op['template']['guided_segmentation']:             segmentation_tg_op['phi_propotion_cutoff'] = op['template']['guided_segmentation']['phi_propotion_cutoff']

    # else make segmentation file None
    else:
        segmentation_op = None
        segmentation_tg_op = None

    # Checking existance of data file in output directory, if exist load it.
    data_json_file = os.path.join(out_dir, 'data.json')
    if os.path.exists(data_json_file):
        print 'Loading available data file:' + data_json_file + '     under directory:' + out_dir
        print 'If you want to create initial random alignment again, delete that file.'
        with open(data_json_file) as f:     data_json = json.load(f)

    # else create data file with random initial tranformations/alignment
    else:
        if op['random_initial_alignment']:
            print 'Forcing use of random initial orientations'
            # As no initial angles are available, set random ones
            for rec in data_json:
                rec['loc'] = [0.0, 0.0, 0.0]
                rec['angle'] = [_ for _ in N.random.random(3) * (N.pi * 2)] # N.random.random(3) creates three random numbers between [0,1] and multiplying it with 2*pi gives 3 random angles in 3D space
                #rec['angle'] = (N.random.uniform(-N.pi, N.pi, 3)).tolist() # N.random.random(3) creates three random numbers between [0,1] and multiplying it with 2*pi gives 3 random angles in 3D space

            # Write this in output directory, data_json already contains some information about suntomograms, this will copy the file and add random orientation to the data.
            with open(data_json_file, 'wb') as f:       json.dump(data_json, f, indent=2)
        else:
            print "Neither random initial alignment is forced true in parameters list, nor initial random alignments are available in out_dir from any previous runs"

    # get size of subtomogram, read any subtomogram and get's it shape using numpy
    v = IV.get_mrc( data_json[0]['subtomogram'] )
    size = N.array(v.shape)
    del v

    file_stat_file = os.path.join(out_dir, 'file_stat.json')
    if os.path.isfile(file_stat_file):
        with open(file_stat_file) as f:        file_stat = json.load(f)
        # Read existing file_stat file and retrieve passes
        file_stat['passes'] = {int(_):file_stat['passes'][_] for _ in file_stat['passes']}

    else:
        # else create file_stat file which contains information about data file, stat file, out_dir, passes completed and current pass
        # This is used to start the MPP from current pass again, if it stops/terminates in middle due to some issue.
        file_stat = {}
        file_stat['out_dir'] = out_dir
        file_stat['file_stat_file'] = file_stat_file
        file_stat['data_json_file'] = data_json_file     # this is the initial data_json_file
        file_stat['pass_i_current'] = 0
        file_stat['passes'] = {}
        with open(file_stat_file, 'w') as f:        json.dump(file_stat, f, indent=2)


    # Clustering method used in each pass/iteration of pattern mining
    cluster_modes = ['kmeans']

    # Use this to collect all cluster averages and their FSC information for cluster selection purpose. IMPORTANT: we do not store those information that are subject to change to here, for example 'is_specific'
    cluster_info = defaultdict(dict)

    # in order to save storage, we store those information that may subject to change here, for example 'is_specific'
    cluster_info_stat = defaultdict(dict)

    cas_re = None
    fsc_stat = {}

    for pass_i in range(op['options']['pass_num']):

        # we only need data from the last completed pass
        if pass_i < file_stat['pass_i_current']:   continue
        
        # n_chunk is used to control the number of tasks given to single worker
        # either "min chunk size" in parameters list or (no. of tasks / min worker number)... min worker number is also set in parameters list under "options"
        print "Worker Number ", self.runner.work_queue.get_worker_number()
        print "Min chunk size ",  op['options']['min_chunk_size']
        print "Data_length", len(data_json), "Other option is ",  int( N.ceil(len(data_json) / max(op['options']['min_worker_num'], self.runner.work_queue.get_worker_number()+10)) )
        n_chunk = max(op['options']['min_chunk_size'], int( N.ceil(len(data_json) / max(op['options']['min_worker_num'], self.runner.work_queue.get_worker_number()+10)) ))
        print "n_chunk", n_chunk
        
        print '-------------------------------------------------------'
        print 'Beginning pass #%d' % (pass_i)       ;       sys.stdout.flush()
        pass_start_time = time.time()

        # Create pass directory if not present
        pass_dir  = os.path.join(out_dir, 'pass', '%03d' % (pass_i))
        if not os.path.exists(pass_dir) :            os.makedirs(pass_dir)

        # Create files
        data_json_file =  os.path.join(pass_dir, 'data.json')
        cluster_info_file = os.path.join(pass_dir, 'cluster_info.pickle')
        cluster_info_stat_file = os.path.join(pass_dir, 'cluster_info_stat.pickle')
        cluster_modes_file = os.path.join(pass_dir, 'cluster_modes.json') # Store clustering method
        cluster_average_select_file = os.path.join(pass_dir, 'cluster_average_select.pickle')
        align_template_file = os.path.join(pass_dir, 'align_template.pickle')
        fsc_stat_file = os.path.join(pass_dir, 'fsc_stat.json')

        # Say we are at pass number i, now we skipped all passes from 1 to i-1, we want to check whether ith pass was completed or not.
        # if fsc_stat_file in ith pass directory is available that means ith pass was completed and we need to load data from ith pass directoy and continue to (i+1)th pass
        if os.path.isfile(fsc_stat_file):

            # first load fsc_stat to determine whether need to stop
            print 'loading', fsc_stat_file
            with open(fsc_stat_file, 'r') as f:    fsc_stat = json.load(f)
            fsc_stat = fsc_stat_json_convert(fsc_stat) #function call

            if os.path.isfile(cluster_modes_file):
                print 'loading', cluster_modes_file
                with open(cluster_modes_file) as f:     cluster_modes = json.load(f)
            # Because we have skipped all the passes before the latest one where the last run was stopped/terminated
            # We want to see at what step it was stopped, whether it satisfied the actual termination condition or if not what kind of step it was, normal kmeans or adaptive kmeans or sequential expansion, etc.
            # Stop function call in next line takes care of this checking
            sts_re = stop_test_stat(op=op, fsc_stat=fsc_stat, pass_i=pass_i, cluster_modes=cluster_modes)
            # Updating cluster modes, if stop function appended adaptive_kmeans or sequential expansion then it will become that
            cluster_modes = sts_re['cluster_modes']

            if sts_re['should_stop']:
                print 'no more improvements seen, stop'
                break

            # load other files

            print 'loading', data_json_file
            with open(data_json_file) as f:     data_json = json.load(f)

            print 'loading ', cluster_average_select_file
            with open(cluster_average_select_file, 'rb') as f:    cas_re = pickle.load(f)

            print 'loading', cluster_info_stat_file 
            with open(cluster_info_stat_file, 'rb') as f:                    cluster_info_stat = pickle.load(f)

            # load previous clus_info generated from each iteration
            print 'loading cluster info' 
 
            for pass_i_t in file_stat['passes']:
                print 'loading', file_stat['passes'][pass_i_t]['cluster_info_file']
                with open(file_stat['passes'][pass_i_t]['cluster_info_file'], 'rb') as f:                    cluster_info[pass_i_t] = pickle.load(f)

            if os.path.isfile(align_template_file):
                print 'loading ' + align_template_file
                with open(align_template_file, 'rb') as f:    at_ress = pickle.load(f)

                at_ress = [_.result for _ in at_ress]


            print '....  and go to next pass'

            continue



        # If fsc_stat_file is not avialable that means current pass was not completed and we need to start it again.
        assert      pass_i not in file_stat['passes']
        # Because in file_stat we save link to each file in pass folder, we already have data for passes < pass_i, so we start creating new data for current pass_i
        file_stat['passes'][pass_i] = {}
        file_stat['passes'][pass_i]['pass_i'] = pass_i
        file_stat['passes'][pass_i]['pass_dir'] = pass_dir
        # Set pass_i_current to pass_i
        file_stat['pass_i_current'] = pass_i

        # keep the copy of configuration/parameter file
        op_file_t = os.path.join(pass_dir, 'pursuit-op-%d.json'%(int(time.time()), ))
        with open(op_file_t, 'w') as f:        json.dump(op, f, indent=2)
        

        file_stat['passes'][pass_i]['cluster_modes_file'] = cluster_modes_file
        # processing for cluster_modes_file, if it doesn't exist then save it with information in variable "cluster_modes" which represents type of clustering used.
        if os.path.isfile(cluster_modes_file):
            print 'loading', cluster_modes_file
            with open(cluster_modes_file) as f:     cluster_modes = json.load(f)
        else:
            with open(cluster_modes_file, 'w') as f:     json.dump(cluster_modes, f)


        # dictionary of subtomograms in data_json, just subtomogram file links, nothing else
        #data_json_copy = copy.deepcopy(data_json)
        if pass_i == 0:
            cluster_dump_file = os.path.join(pass_dir, 'cluster.pickle')
            file_stat['passes'][pass_i]['cluster_dump_file'] = cluster_dump_file
            # if clusters.pickle available then load it else we will create it by doind dimension reduction and kmeans clustering
            if os.path.exists(cluster_dump_file):
                print 'loading', cluster_dump_file
                with open(cluster_dump_file, 'rb') as f:
                    cluster_dump_file__load = pickle.load(f)
                fsc_dict = cluster_dump_file__load['fsc_dict']
                clusters_populated = cluster_dump_file__load['clusters_populated']
            
            else:
                cluster_kmeans_dump_file = os.path.join(pass_dir, 'cluster__kmeans.pickle')
                file_stat['passes'][pass_i]['cluster_kmeans_dump_file'] = cluster_kmeans_dump_file
                clusters_populated = {}              # clusters_populated for current pass
                fsc_dict = {}           # collect fsc
                # If kmeans file doesn't exist then create one
                start_time = time.time()
                if not os.path.isfile(cluster_kmeans_dump_file):    
                    kmeans_clusters = {}
                    kmeans_labels = N.array((), dtype=N.int)
                    for pass_0 in range(op['options']['num_of_first_iterations']):
                        #randomize angles in data_json
                        data_json_copy = copy.deepcopy(data_json)
                        if pass_0 == 0:
                            temp_seed = N.random.random(1)
                        N.random.seed(int(100000*temp_seed))
                        for rec,dj in enumerate(data_json_copy):
                            dj['angle'] = [_ for _ in N.random.random(3) * (N.pi * 2)]
                            temp_seed = N.random.random(1)
                                                                            
                        #print data_json_copy[0]['angle'], data_json_copy[0]['subtomogram'][(len(data_json_copy[0]['subtomogram'])-10):len(data_json_copy[0]['subtomogram'])]
                        data_json_dict = {_['subtomogram']:_ for _ in data_json_copy}
                        temp_file_name = 'dimension_reduction_' + str(pass_0) + '.pickle'
                        dimension_reduction_dump_file = os.path.join(pass_dir, temp_file_name)
                        file_stat['passes'][pass_i]['dimension_reduction_dump_file_'+str(pass_0)] = dimension_reduction_dump_file
                        # if dim red file exist load it else do dim red
                        if os.path.exists(dimension_reduction_dump_file):
                            print 'loading ' + dimension_reduction_dump_file
                            with open(dimension_reduction_dump_file, 'rb') as f:
                                dimension_reduction_dump_file__load = pickle.load(f)
                            cfp_re = dimension_reduction_dump_file__load['cfp_re']

                        else:
                            print 'Dimension reduction for subpass-' + str(pass_0)
                            dim_start_time = time.time()
                            data_json_pca_train = None
                            if (op['dim_reduction']['train_with_selected_clusters_only']) and (cas_re is not None):
                                data_json_pca_train__set = set()
                                data_json_pca_train = []
                                for c in cas_re['selected_templates']:
                                    t = cas_re['tk_info'][cas_re['selected_templates'][c]['subtomogram']]
                                    if op['dim_reduction']['restrict_to_specific_clusters'] and (cluster_info_stat[t['pass_i']][t['cluster']]['is_specific'] is not None):
                                        continue #continue if it is non-specific
                                    for d in t['data_json']:
                                        if d['subtomogram'] in data_json_pca_train__set:
                                            continue
                                        data_json_pca_train__set.add(d['subtomogram'])
                                        data_json_pca_train.append(data_json_dict[d['subtomogram']])
                                del c,t
                            if op['dim_reduction']['with_missing_wedge']:
                                cfp_re = CU.covariance_filtered_pca_with_wedge(self=self, data_json=data_json_copy, n_chunk=n_chunk, op=op, pass_dir=pass_dir, max_feature_num=op['dim_reduction']['max_feature_num'])
                            else:
                                cfp_re = CU.covariance_filtered_pca(self=self, data_json_model=data_json_pca_train, data_json_embed=data_json_copy, normalize=op['dim_reduction']['normalize'], segmentation_tg_op=(segmentation_tg_op if op['dim_reduction']['use_segmentation_mask'] else None), n_chunk=n_chunk, pca_op=op['dim_reduction']['pca'], max_feature_num=op['dim_reduction']['max_feature_num'])

                            with open(dimension_reduction_dump_file, 'wb') as f:    pickle.dump({'cfp_re':cfp_re, 'data_json_pca_train':data_json_pca_train, 'data_json':data_json_copy}, f, protocol=-1)
                            print "Dimension Reduction took: %2.6f sec" % (time.time() - dim_start_time)

                        sys.stdout.flush()
                        # Dimesion reduction done.
                        print "Generating subtomogram sets through kmeans: subpass-" + str(pass_0)
                        kmeans_k = op['cluster']['kmeans']['number']
                        print 'k =', kmeans_k
                        kmeans_labels_temp = CU.kmeans_clustering(x=cfp_re['red'], k=kmeans_k)
                        kmeans_labels_temp = kmeans_labels_temp + (pass_0 * kmeans_k)
                        kmeans_labels = N.concatenate((kmeans_labels,kmeans_labels_temp),0)
                        kmeans_clusters_temp = CU.labels_to_clusters(data_json=data_json_copy, labels=kmeans_labels_temp, cluster_mode=('kmeans' if 'kmeans-adaptave' not in set(cluster_modes) else 'kmeans-adaptave'))
                        kmeans_clusters.update(kmeans_clusters_temp)
                    
                    # calculate ssnr
                    cluster_ssnr_fsc__op = {}
                    cluster_ssnr_fsc__op['ssnr'] = copy.deepcopy(op['ssnr'])
                    if op['cluster']['ssnr']['segmentation']:
                        cluster_ssnr_fsc__op['segmentation_tg'] = copy.deepcopy(segmentation_tg_op)
                    cluster_ssnr_fsc = CU.cluster_ssnr_fsc(self=self, clusters={_:kmeans_clusters[_]['data_json'] for _ in kmeans_clusters}, n_chunk=n_chunk, op=cluster_ssnr_fsc__op)
                    with open(cluster_kmeans_dump_file, 'wb') as f:
                        pickle.dump({'cluster_ssnr_fsc':cluster_ssnr_fsc, 'kmeans_k':kmeans_k*op['options']['num_of_first_iterations'],'kmeans_labels':kmeans_labels, 'kmeans_clusters':kmeans_clusters}, f, protocol=-1)

                # if kmeans file exist then load it
                else:
                    print 'loading', cluster_kmeans_dump_file
                    with open(cluster_kmeans_dump_file, 'rb') as f:     tmp = pickle.load(f)
                    cluster_ssnr_fsc = tmp['cluster_ssnr_fsc']
                    kmeans_labels = tmp['kmeans_labels']
                    kmeans_clusters = tmp['kmeans_clusters']
                    del tmp

                label_t = N.max([_ for _ in clusters_populated]) + 1 if (len(clusters_populated) > 0) else 0
                for c in kmeans_clusters:
                    clusters_populated[label_t] = kmeans_clusters[c]
                    clusters_populated[label_t]['original_label'] = c
                    fsc_dict[label_t] = cluster_ssnr_fsc['fsc'][c]
                    label_t += 1

                assert  len(clusters_populated) > 0
                with open(cluster_dump_file, 'wb') as f:    pickle.dump({'fsc_dict':fsc_dict, 'clusters_populated':clusters_populated}, f, protocol=-1)
                print 'Number of generated averages', len(clusters_populated)
                print "Subtomogram cluster generation and FSC calculation  : %2.6f sec" % (time.time() - start_time)
        
        else: #if pass_i is not 0
            #data_json = copy.deepcopy(data_json_backup)
	    data_json_dict = {_['subtomogram']:_ for _ in data_json}

            sys.stdout.flush()

            #--------------------
	    # Dimension Reduction Step
	    # create dimension reduction file
	    dimension_reduction_dump_file = os.path.join(pass_dir, 'dimension_reduction.pickle')
	    # write in file_stat
	    file_stat['passes'][pass_i]['dimension_reduction_dump_file'] = dimension_reduction_dump_file

	    # IF dimension reduction file exist load it else create one.
	    if os.path.exists(dimension_reduction_dump_file):
	        print 'loading ' + dimension_reduction_dump_file
    		with open(dimension_reduction_dump_file, 'rb') as f:    dimension_reduction_dump_file__load = pickle.load(f)
		cfp_re = dimension_reduction_dump_file__load['cfp_re']

	    # creating dimension reduction file
	    else:
                dim_start_time = time.time()
		data_json_pca_train = None
		# For first iteration car_re will be None because we don't have any selected clusters from last iteration
		# for further iterations it contains information about subtomograms that were part of selected patterns, so that only these subtomograms can be used to reduce dimensions
		if (op['dim_reduction']['train_with_selected_clusters_only']) and (cas_re is not None):
		    # use the subtomograms of the selected clusters from last iteration (with latest alignment) for dimension reduction'
		    #'only specific clusters' if op['dim_reduction']['restrict_to_specific_clusters']
		        
		    data_json_pca_train__set = set()
		    data_json_pca_train = []
		    for c in cas_re['selected_templates']:
		        t = cas_re['tk_info'][cas_re['selected_templates'][c]['subtomogram']] # this select subtomogram and it's further information about it stored in tk_info

		        # Selecting specific clusters if "restrict to specific cluster" is true. This information is selected from last iteration
		        if op['dim_reduction']['restrict_to_specific_clusters'] and (cluster_info_stat[t['pass_i']][t['cluster']]['is_specific'] is not None):  continue #continue if it is non-specific

		        # If restrict to specific clusters is false, if it was specific
		        for d in t['data_json']: # data_json is key inside subtomogram in tk_info, this list is of subtomograms belonging to pattern "c"
		            if d['subtomogram'] in data_json_pca_train__set:    continue # continue if subtomogram was already in training set, due to overlap in patterns
		            # Other wise add to set. This adds all the info about subtomogram
		            data_json_pca_train__set.add(d['subtomogram'])
		            # This just adds the name of subtomogram from list data_json_dict
		            data_json_pca_train.append(data_json_dict[d['subtomogram']]) 
		            # Add the corresponding subtomogram of specific clusters selected from the last iteration, with latest alignment information.
		            # The use of latest alignment info is because we want to use dimension reduction and clustering to obtain new and better clusters

		    del c, t

		# Pass all the information to dimension reduction function depending on missing wedge parameter
		# Calculates dimensions by first imputing all the subtomograms, then calculating covariance and cov avg, then using cutoff for top 10000 covariances and calculating desired number of dimensions using EM-PCA
		if op['dim_reduction']['with_missing_wedge']:
		    cfp_re = CU.covariance_filtered_pca_with_wedge(self=self, data_json=data_json, n_chunk=n_chunk, op=op, pass_dir=pass_dir, max_feature_num=op['dim_reduction']['max_feature_num'])
		else:
		    cfp_re = CU.covariance_filtered_pca(self=self, data_json_model=data_json_pca_train, data_json_embed=data_json, normalize=op['dim_reduction']['normalize'], segmentation_tg_op=(segmentation_tg_op if op['dim_reduction']['use_segmentation_mask'] else None), n_chunk=n_chunk, pca_op=op['dim_reduction']['pca'], max_feature_num=op['dim_reduction']['max_feature_num'])

		# cfp_re contains dataset names 'red' which is num_of_subtomograms * 50 table, i.e. value of coordinates of all subtomorgams in 50 coordinate system deducted using PCA
		# Save all the information about dimension redunction step, reduced dimensions returned from function, subtomograms used for reduction and whole data_json
		with open(dimension_reduction_dump_file, 'wb') as f:    pickle.dump({'cfp_re':cfp_re, 'data_json_pca_train':data_json_pca_train, 'data_json':data_json}, f, protocol=-1)

		print "Dimension Reduction took: %2.6f sec" % (time.time() - dim_start_time)

	    sys.stdout.flush()

	    # -------------------------------------------
	    # generate clusters

	    # Creating another file, could be created earlier with other files as well.
	    cluster_dump_file = os.path.join(pass_dir, 'cluster.pickle')
	    # Setting link in file_stat file
	    file_stat['passes'][pass_i]['cluster_dump_file'] = cluster_dump_file

	    # If cluster.pickle already exist then load it else create one
	    if os.path.exists(cluster_dump_file):
	        print 'loading', cluster_dump_file
	        
	        with open(cluster_dump_file, 'rb') as f:        cluster_dump_file__load = pickle.load(f)
	        fsc_dict = cluster_dump_file__load['fsc_dict']
	        clusters_populated = cluster_dump_file__load['clusters_populated']          # clusters_populated for current pass

	    else:
	        start_time = time.time()

	        clusters_populated = {}              # clusters_populated for current pass
	        fsc_dict = {}           # collect fsc

		if 'kmeans' in set(cluster_modes):            
		    #   now we generate subtomogram sets through kmeans
		    cluster_kmeans_dump_file = os.path.join(pass_dir, 'cluster__kmeans.pickle')
		    file_stat['passes'][pass_i]['cluster_kmeans_dump_file'] = cluster_kmeans_dump_file

		    # File doesn't exist then create else load directly.
		    if not os.path.isfile(cluster_kmeans_dump_file):
		        print   'Generate subtomogram sets through kmeans clustering',
		        kmeans_k = op['cluster']['kmeans']['number']

		        # Check whether we need to do kmeans-adaptive now or just kmeans
		        # if cluster_modes include kmeans-adaptive then we instead of using direct kmeans number = 10, we use kmeans-adaptive-ratio * specific cluster count
		        if 'kmeans-adaptave' not in set(cluster_modes):
		            kmeans_k = op['cluster']['kmeans']['number']

		        else:

		            specific_cluster_count = 0

		            for c in cas_re['selected_templates']:
		                t = cas_re['tk_info'][cas_re['selected_templates'][c]['subtomogram']]
		                if cluster_info_stat[t['pass_i']][t['cluster']]['is_specific'] is not None:  continue        # if this cluster is selected in last iteration, but identified as non-specific, ignore it
		                specific_cluster_count += 1
	     
		            kmeans_k = int(N.round(specific_cluster_count * float(op['cluster']['kmeans']['adaptive_k_ratio'])))
		            del c, t, specific_cluster_count

		        print 'k =', kmeans_k

		        # function call below returns labels (representing cluster) for each data point i.e. subtomogram
		        kmeans_labels = CU.kmeans_clustering(x=cfp_re['red'], k=kmeans_k)
		        # divides data according to labels
		        kmeans_clusters = CU.labels_to_clusters(data_json=data_json, labels=kmeans_labels, cluster_mode=('kmeans' if 'kmeans-adaptave' not in set(cluster_modes) else 'kmeans-adaptave'))
		        # Calculate ssnr
		        # ssnr for each cluster is array size (x/2 + 1), where x is miminum among subtomogram dimensions
		        cluster_ssnr_fsc__op = {}
		        cluster_ssnr_fsc__op['ssnr'] = copy.deepcopy(op['ssnr'])
		        if op['cluster']['ssnr']['segmentation']:       cluster_ssnr_fsc__op['segmentation_tg'] = copy.deepcopy(segmentation_tg_op)
		        cluster_ssnr_fsc = CU.cluster_ssnr_fsc(self=self, clusters={_:kmeans_clusters[_]['data_json'] for _ in kmeans_clusters}, n_chunk=n_chunk, op=cluster_ssnr_fsc__op)
		        # Write this in cluster_kmeans.pickle file
		        with open(cluster_kmeans_dump_file, 'wb') as f:    pickle.dump({'cluster_ssnr_fsc':cluster_ssnr_fsc, 'kmeans_k': kmeans_k,'kmeans_labels':kmeans_labels, 'kmeans_clusters':kmeans_clusters}, f, protocol=-1)

		    else:
		        print 'loading', cluster_kmeans_dump_file

		        with open(cluster_kmeans_dump_file, 'rb') as f:     tmp = pickle.load(f)
		        cluster_ssnr_fsc = tmp['cluster_ssnr_fsc']
		        kmeans_labels = tmp['kmeans_labels']
		        kmeans_clusters = tmp['kmeans_clusters']
		        del tmp

		    label_t = N.max([_ for _ in clusters_populated]) + 1 if (len(clusters_populated) > 0) else 0
		    for c in kmeans_clusters:
		        clusters_populated[label_t] = kmeans_clusters[c]        ;       clusters_populated[label_t]['original_label'] = c
		        fsc_dict[label_t] = cluster_ssnr_fsc['fsc'][c]
		        label_t += 1
		# Finally clusters_populated contains info about subtomograms in each cluster
		# fsc_stat contains info about fsc scores of each cluster
		    
		     
		# Sequential expansion happens after few iterations when kmeans can't give any improvements
		if 'sequential' in set(cluster_modes):
		    # the sequential expansion is independent to dimension reduction. It comes directly from the alignment and segmentation of last iteration
		    cluster_sequential_dump_file = os.path.join(pass_dir, 'cluster__sequential.pickle')
		    file_stat['passes'][pass_i]['cluster_sequential_dump_file'] = cluster_sequential_dump_file

		    if not os.path.isfile(cluster_sequential_dump_file):
		        print 'Generate subtomogram sets through sequential expansion'

		        # at_ress is alignment file
		        al_t = [_['align'] for _ in at_ress]
		        if ('filtering__second_largest_cut' in op['cluster']['sequential_expansion']) and (len(al_t[0]) > 1):
		            # when there are multiple templates
		            import tomominer.pursuit.multi.recursive.filtering.second_largest_cut as PMRFS
		            data_json_sp = PMRFS.do_filter(al=al_t, dj=data_json)
		        else:
		            # when there is only a single template
		            data_json_sp = data_json
		     
		        pmcra_op = copy.deepcopy(op['cluster']['sequential_expansion'])
		        pmcra_op['min_expansion_size'] = op['cluster']['size_min']
		        pmcra_op['n_chunk'] = n_chunk
		        pmcra_op['ssnr_sequential'] = {}
		        pmcra_op['ssnr_sequential']['ssnr'] = copy.deepcopy(op['ssnr'])
		        if op['cluster']['ssnr']['segmentation']:   pmcra_op['ssnr_sequential']['segmentation_tg'] = copy.deepcopy(segmentation_tg_op)
		        pmcra_re = CU.cluster_formation_alignment_fsc__by_global_maximum(self=self, dj=data_json_sp, op=pmcra_op)
		        pmcra_re['data_json'] = data_json_sp
		          
		        with open(cluster_sequential_dump_file, 'wb') as f:    pickle.dump(pmcra_re, f, protocol=-1)

		    else:
		        # this step sometimes is very time consuming, so we just load exising file

		        print 'loading', cluster_sequential_dump_file

		        with open(cluster_sequential_dump_file, 'rb') as f:    pmcra_re = pickle.load(f)

		    label_t = N.max([_ for _ in clusters_populated]) + 1 if (len(clusters_populated) > 0) else 0
		    for c in pmcra_re['dj_gm']:                    
		        clusters_populated[label_t] = {'cluster_mode':'sequential', 'data_json':pmcra_re['dj_gm'][c]['data_json'], 'original_label':c}
		        fsc_dict[label_t] = pmcra_re['dj_gm'][c]['fsc']
		        label_t += 1



		assert  len(clusters_populated) > 0
		# Saving cluster.pickle file that contains information about clusters (subtomograms in each cluster and fsc score)
		with open(cluster_dump_file, 'wb') as f:    pickle.dump({'fsc_dict':fsc_dict, 'clusters_populated':clusters_populated}, f, protocol=-1)


		print 'Number of generated averages', len(clusters_populated)
		print "Subtomogram cluster generation and FSC calculation  : %2.6f sec" % (time.time() - start_time)

        #----------------------------------------- 
        # remove small clusters for further processing
        # although all the clusters small or large, are saved in cluster.pickle file for this pass
        for c in list(clusters_populated.iterkeys()):
            if len(clusters_populated[c]['data_json']) < op['cluster']['size_min']:      del clusters_populated[c]

        #----------------------------------
        # update cluster_info and cluster_info_stat file 
        # cluster_info file contains information about each cluster like: pass_i, template, fsc, data_json, cluster_mode
        # Whereas cluster_info_stat contains info about each pass and for each pass it contains info about cluster size and specificity for each cluster
        for c in clusters_populated:
            if c not in cluster_info[pass_i]:       cluster_info[pass_i][c] = {}
            cluster_info[pass_i][c]['pass_i'] = pass_i
            cluster_info[pass_i][c]['cluster'] = c
            if c in fsc_dict:   cluster_info[pass_i][c]['fsc'] = fsc_dict[c]
            cluster_info[pass_i][c]['cluster_mode'] = clusters_populated[c]['cluster_mode']
            cluster_info[pass_i][c]['data_json'] = clusters_populated[c]['data_json']


        for c in cluster_info[pass_i]:
            if pass_i not in cluster_info_stat:     cluster_info_stat[pass_i] = {}
            if c not in cluster_info_stat[pass_i]:       cluster_info_stat[pass_i][c] = {'cluster_size':len(cluster_info[pass_i][c]['data_json'])}


        sys.stdout.flush()


        #--------------------------------------------------------------------
        # Using the clusters we just generated, we build cluster averages.

        cluster_averaging_file = os.path.join(pass_dir, 'cluster_averaging.pickle')
        file_stat['passes'][pass_i]['cluster_averaging_file'] = cluster_averaging_file

        # if cluster average file exist load it else create
        if os.path.exists(cluster_averaging_file):
            print 'loading ' + cluster_averaging_file
            with open(cluster_averaging_file, 'rb') as f:    ca_re = pickle.load(f)
        else:
            print 'Calculate subtomogram averages'
            cu_ca_op = copy.deepcopy(op['cluster']['averaging'])
            cu_ca_op['out_dir'] = os.path.join(pass_dir, 'clus_avg')
            cu_ca_op['pass_i'] = pass_i         # this is used to be recorded into the template_key
            cu_ca_op['n_chunk'] = n_chunk
            if 'smooth' in cu_ca_op:                cu_ca_op['smooth']['fsc'] = {_:fsc_dict[_] for _ in fsc_dict}
            # function call for calculating cluster averages
            # ca_re will contain template_keys as in cluster_averaging file
            ca_re = CU.cluster_averaging(self=self, clusters={_:clusters_populated[_]['data_json'] for _ in clusters_populated}, op=cu_ca_op)
            with open(cluster_averaging_file, 'wb') as f:    pickle.dump(ca_re, f, protocol=1)

        template_keys = ca_re['template_keys']
        #print template_keys

        # update cluster_info with template_keys
        for c in template_keys:
            if c not in cluster_info[pass_i]:       cluster_info[pass_i][c] = {}
            cluster_info[pass_i][c]['template_key'] = template_keys[c]
            assert  cluster_info[pass_i][c]['template_key']['cluster'] == c
            assert  cluster_info[pass_i][c]['template_key']['pass_i'] == pass_i


        sys.stdout.flush()


        #-------------------------------------------------------------------
        # select clusters
        # Load select cluster average file if available else create it

        file_stat['passes'][pass_i]['cluster_average_select_file'] = cluster_average_select_file
        if os.path.exists(cluster_average_select_file):
            print 'loading ' + cluster_average_select_file
            with open(cluster_average_select_file, 'rb') as f:    cas_re = pickle.load(f)
        else:
            print 'select averages'
            select_op = copy.deepcopy(op['cluster_average_select_fsc'])
            select_op['cluster'] = op['cluster']
            select_op['align_op'] = op['align']
            cas_re = CU.cluster_average_select_fsc(self=self, cluster_info=cluster_info, cluster_info_stat=cluster_info_stat, op=select_op)

            if op['common_frame']['mode'] == 1:
                caacfmp = CU.cluster_average_align_common_frame__multi_pair(self=self, tk=cas_re['selected_templates'], align_op=op['align'], loc_r_max=size.min()*float(op['common_frame']['loc_r_max_proportion']), pass_dir=pass_dir)
                cas_re['selected_templates_common_frame'] = caacfmp['tka']
                cas_re['selected_templates_common_frame__unrotated_clus'] = caacfmp['unrotated_clus']
                cas_re['selected_templates_common_frame__align_to_clus'] = caacfmp['align_to_clus']
                cas_re['cluster_average_align_common_frame__multi_pair_file'] = os.path.join(pass_dir, 'cluster_average_align_common_frame__multi_pair.pickle')
                with open(cas_re['cluster_average_align_common_frame__multi_pair_file'], 'wb') as f:    pickle.dump(caacfmp, f, protocol=-1)

            elif op['common_frame']['mode'] == 2:
                caacfmp = CU.cluster_average_align_common_frame__single_best(self=self, tk=cas_re['selected_templates'], align_op=op['align'], pass_dir=pass_dir)
                cas_re['selected_templates_common_frame'] = caacfmp['tka']
                cas_re['cluster_average_align_common_frame__single_best_file'] = os.path.join(pass_dir, 'cluster_average_align_common_frame__single_best.pickle')
                with open(cas_re['cluster_average_align_common_frame__single_best_file'], 'wb') as f:    pickle.dump(caacfmp, f, protocol=-1)


            else:
                raise   Exception('op[\'common_frame\'][\'mode\']'%(op['common_frame']['mode'],))

            # segment those templates aligned to the common frame and store the segmented files. No information is aved about them 
            if segmentation_op is not None:
                print 'segmenting the selected and aligned averages'
                template_segmentation_op = copy.deepcopy(segmentation_op)
                if ('segmentation' in op['template']) and ('normalize_and_take_abs' in op['template']['segmentation']) and (op['template']['segmentation']['normalize_and_take_abs']):
                    template_segmentation_op['normalize_and_take_abs'] = True
                CU.template_segmentation(self=self, tk=cas_re['selected_templates_common_frame'], op=template_segmentation_op)

            with open(cluster_average_select_file, 'wb') as f:    pickle.dump(cas_re, f, protocol=-1)

        if True:
            # print information on selected templates
            common_path_prefix = os.path.commonprefix([  cas_re['selected_templates'][_]['subtomogram'] for _ in cas_re['selected_templates']  ])
            print len(cas_re['selected_templates']), 'averages selected.'
            print 'Average list with common path prefix:', common_path_prefix
            print 'average id', '\t', 'generation mode', '\t', 'SFSC score', '\t', 'set size', '\t', 'average file'
            for _ in cas_re['selected_templates']:       
                print _, '\t\t', cas_re['tk_info'][cas_re['selected_templates'][_]['subtomogram']]['cluster_mode'], '\t\t', cas_re['tk_info'][cas_re['selected_templates'][_]['subtomogram']]['fsc'].sum(), '\t', len(cas_re['tk_info'][cas_re['selected_templates'][_]['subtomogram']]['data_json']), '\t', cas_re['selected_templates'][_]['subtomogram'][len(common_path_prefix):]

        #sys.stdout.flush()


 
        #----------------------------------------------------------
        # align against selected templates
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

            at_ress = CU.align_to_templates__batch(self=self, op=op, data_json=data_json, segmentation_tg_op=segmentation_tg_op, tmp_dir=align_template__tmp_dir, tem_keys=cas_re['selected_templates_common_frame'])

            with open(align_template_file, 'wb') as f:    pickle.dump(at_ress, f, protocol=-1)
            shutil.rmtree(align_template__tmp_dir)          #  clean up temporary alignment info

            print "Align all volumes to cluster_templates. %2.6f sec" % (time.time() - start_time)

        at_ress = [_.result for _ in at_ress]


        sys.stdout.flush()



        #--------------------------------------------------------------------------
        # select templates according to matching specificity, and re-calculate best alignments only according to the selected templates
        cratcms = CU.cluster_removal_according_to_center_matching_specificity(ci=cluster_info, cis=cluster_info_stat, al=at_ress, tk=cas_re['selected_templates'], significance_level=op['cluster_removal_according_to_center_matching_specificity']['significance_level'])
        with open(os.path.join(pass_dir, 'cluster_removal_according_to_center_matching_specificity.pickle'), 'wb') as f:     pickle.dump(cratcms, f, protocol=-1)


        with open(cluster_info_file, 'wb') as f:    pickle.dump(cluster_info[pass_i], f, protocol=-1)           # we only dump cluster info of CURRENT PASS, to save storage!!
        file_stat['passes'][pass_i]['cluster_info_file'] = cluster_info_file

        with open(cluster_info_stat_file, 'wb') as f:    pickle.dump(cluster_info_stat, f, protocol=-1)
        file_stat['passes'][pass_i]['cluster_info_stat_file'] = cluster_info_stat_file
        

        # record template alignment results
        data_json_new = []
        for res in at_ress:

            rec = {}
            rec['subtomogram'] = res['vol_key']['subtomogram']
            rec['mask'] = res['vol_key']['mask']

            # numpy arrays are not pickle-able.  Convert to a list of numbers.
            rec['angle'] = [_ for _ in res['best']['angle']]
            rec['loc'] = [_ for _ in res['best']['loc']]
            rec['score'] = res['best']['score']

            if res['best']['template_id'] is not None:
                rec['template'] = cas_re['selected_templates_common_frame'][res['best']['template_id']]     # store the matched template for interpolation of next iteration. Also since cas_re['selected_templates_common_frame'] may contain segmentation information, store a copy to data_json to facilate segmentation of the next iteration
                assert  'id' in rec['template']

            
            data_json_new.append(rec)

        data_json = data_json_new
        with open(data_json_file, 'w') as f:    json.dump(data_json, f, indent=2)
        file_stat['passes'][pass_i]['data_json_file'] = data_json_file

        #print 'Failure best alignment number', len([_ for _ in data_json if not N.isfinite(_['score'])])


        # ----------------------------------------------------
        # calculate fsc_stat for determining stopping criterion
        file_stat['passes'][pass_i]['fsc_stat_file'] = fsc_stat_file
        if os.path.exists(fsc_stat_file):
            print 'loading' + fsc_stat_file
            with open(fsc_stat_file) as f:    fsc_stat = json.load(f)
            fsc_stat = fsc_stat_json_convert(fsc_stat)

        else:

            # calculate FSC sum score of clusters selected in current pass
            fsc_stat_t = {}
            for _ in cas_re['selected_templates']:
                kt = cas_re['selected_templates'][_]['subtomogram']
                fsc_stat_t[_] = { 'pass_i': cas_re['tk_info'][kt]['pass_i'],   'cluster': cas_re['tk_info'][kt]['cluster'],    'fsc': cas_re['tk_info'][kt]['fsc'].sum()          }
                fsc_stat_t[_]['cluster_mode'] = cluster_info[ fsc_stat_t[_]['pass_i'] ][ fsc_stat_t[_]['cluster'] ]['cluster_mode']
            fsc_stat[pass_i] = fsc_stat_t
            del fsc_stat_t, kt, _

            # record specificity of all clusters
            for pass_i_t, fsc_stat_t in fsc_stat.iteritems():
                for template_i, fsc_stat_tt in fsc_stat_t.iteritems():
                    fsc_stat_tt['is_specific'] = cluster_info_stat[fsc_stat_tt['pass_i']][fsc_stat_tt['cluster']]['is_specific']
            del pass_i_t, fsc_stat_t, template_i, fsc_stat_tt

            with open(fsc_stat_file, 'w') as f:    json.dump(fsc_stat, f, indent=2)

        with open(file_stat_file, 'w') as f:      json.dump(file_stat, f, indent=2)


        #--------------------------------------------------------
        # stopping criteron test
        sts_re = stop_test_stat(op=op, fsc_stat=fsc_stat, pass_i=pass_i, cluster_modes=cluster_modes)
        cluster_modes = sts_re['cluster_modes']


        if sts_re['should_stop']:
            print 'no more improvements seen, stop'
            break

        print "Entire pass: %2.6f sec" % (time.time() - pass_start_time)

        sys.stdout.flush()


    return file_stat




'''

todo: change terms: cluster --> pattern,    template --> average

todo: Speedup MPP by removing redundant calculations: If one pattern P has been selected, then we will have all subtomograms aligned with it. In another iteration, if P is selected again, then we DO NOT need to align all subtomograms against it again. This will largely reduce computation time !!! However, the we have to align the subtomograms directly against the original average (also segment if necessary), not the average that is re-aligned into common frames!!!

todo: Also, in a new iteration if the PCA use the same set of subtomograms with same rigid transforms, we can just directly use the PCA result of last iteration.


todo: use UUID to index subtomograms, and use indexes to retrive subtomograms

'''

