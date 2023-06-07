
# iterative refinement using following steps:

# Template guided segmentation and alignment (can add gradient based alignment step. Together with some bandpass mask). Then averaging with missing values.


def average(op, data, out_dir, runner):

    if 'test' in op:
        if ('sample_num' in op['test']) and (len(data) > op['test']['sample_num']):
            print 'testing the procedure using a subsample of %d subtomograms'%(op['test']['sample_num'])
            data = random.sample(data, op['test']['sample_num'])


    # determine the shape of the data.
    v = tomo.parse_mrc(data[0][0])
    vol_shape = v.shape

    avg_key = None
    avg_mask_key = None
    
    ws = None
    old_ws_file = None
    for pass_i in range(op['options']['pass_num']):

        print "Beginning pass #%d" % (pass_i)
        pass_start_time = time.time()

        pass_dir  = os.path.join(out_dir, 'pass_%03d' % (pass_i))


        if not os.path.exists(pass_dir) : 
            os.mkdir(pass_dir);   
            #os.chmod(pass_dir, 0775)    # mxu: enable writing access from other members in the group
            
        ws_file = os.path.join(pass_dir,'dump.pickle')
        if os.path.exists(ws_file):
            print 'jump to next pass'
            old_ws_file = ws_file
            continue
        else:
            if old_ws_file is not None:
                with open(old_ws_file) as f:    ws = pickle.load(f)
                data = ws['data']    


        # Use the selected  cluster to  build cluster average
        gavg_re = avgu.global_average_parallel(runner=runner, data=data, vol_shape=vol_shape, n_chunk=max(op['options']['min_chunk_size'], int(len(selected_subtomograms) / op['options']['worker_num'])), pass_dir=pass_dir, use_fft_avg=op['average']['use_fft_avg'])
        avg_key = gavg_re['vol_avg_out_key']
        avg_mask_key = gavg_re['mask_avg_out_key']

        # see if there is a large bias in missing wedge region
        mask_avg__min = uv.wedge_mask_min(  tomo.parse_mrc(gavg_re['mask_avg_out_key']) )
        print 'mask_avg__min %f'%(mask_avg__min)
        if (pass_i == 0) and (mask_avg__min  <  op['wedge']['avg_min_threshold']):
                print 'stop averaging due to missing wedge'
                return None

        # align all subtomograms against the average
        start_time = time.time()
        at_ress = alu.align_to_templates__parallel(vmal_org=data, vm_tem=[(avg_key, avg_mask_key)], L=op['align']['L'], runner=runner, n_chunk=max(op['options']['min_chunk_size'], int(len(data) / op['options']['worker_num'])), tem_cutoff=0.01)

        new_data = []
        for at_res in at_ress:
            new_data.append( (at_res['vol_key'], at_res['vol_mask_key'], at_res['angle'], at_res['loc']) )



        print 'Align all volumes to cluster_templates. %2.6f sec' % (time.time() - start_time)

        data = new_data



        # record alignment and clustering results           
        data_json = []
        for res in at_ress:
            json_rec = {}
            json_rec['subtomogram'] = res['vol_key']
            json_rec['mask'] = res['vol_mask_key']

            # numpy arrays are not pickle-able.  Convert to a list of numbers.
            json_rec['angle'] = [_ for _ in res['angle']]
            json_rec['loc'] = [_ for _ in res['loc']]
            json_rec['score'] = res['score']
               
            data_json.append(json_rec)

        # write new data to disk
        data_json_file =  os.path.join(pass_dir, 'data_config.json'%(pass_i))
        with open(data_json_file, 'w') as f:
            json.dump(data_json, f, indent=2)

        print "Entire pass: %2.6f sec" % (time.time() - pass_start_time)
        print "---------------------------"

        
        # record and save workspace
        ws = {'pass_i':pass_i, 'pass_i':pass_i, 'data':data, 'gavg_re':gavg_re, 'cfp_re':cfp_re, 'hc_re':hc_re, 'hc_info':hc_info, 'hc_info_f':hc_info_f, 'cs_re':cs_re, 'selected_subtomograms':selected_subtomograms, 'data_json':data_json}

        with open(ws_file, 'w') as f:     pickle.dump(ws, f)



    if ws is None:
        with open(ws_file) as f:    ws = pickle.load(f)

    return ws








if __name__ == '__main__':

    import sys
    qhost = sys.argv[1]
    qport = 5011

    param_file = sys.argv[2]
    data_file  = sys.argv[3]
    out_dir = sys.argv[4]

    with open(param_file) as f:
        op = json.load(f)

    with open(data_file) as f:
        data_json = json.load(f)
    

    data = config.parse_data(data_json)
    
    runner = QueueMaster(qhost, qport)

    average(op=op, data=data, out_dir=out_dir, runner=runner)


