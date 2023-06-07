import json, os, gc

def main():
    with open("merge_all_info_json_files__config.json") as f:    opt = json.load(f)
    with open('model_tomogram_subtomogram__config.json') as f:    op = json.load(f)
    with open(op['model_generation_config_file']) as f:     model_op = json.load(f)
    with open(op['back_projection_reconstruction_config_file']) as f:     bp_op = json.load(f)

    snr = str(opt['SNR'])
    #info_copy_json = opt["file_with_all_subtomograms"]
    #info_json = opt["file_with_all_subtomograms_that_have_pdb_id"]
    
    processing_dir_root = os.path.join(op['out_dir'], str(model_op['copy_number']['total']), '%d-%d-%d'%(model_op['packing']['param']['box']['x'], model_op['packing']['param']['box']['y'], model_op['packing']['param']['box']['z']))
    
    info_copy_json = []
    for model_id in range(model_op["model_number"]):
        if model_id not in op["selected_models"]:    continue
        indir = os.path.abspath( os.path.join(processing_dir_root, str(model_id), str(bp_op['model']['missing_wedge_angle']), snr, ('peak-dog' if 'peak' in op else 'peak-true'), str(op['subtomogram']['size_ratio'] if 'size_ratio' in op['subtomogram'] else op['subtomogram']['size']), 'subtomograms'))
        info_file = os.path.join(indir, "info.json")
        with open(info_file) as f:    info = json.load(f)
        print "#subtomograms in", model_id, "is", len(info)
        info_copy_json.extend(info)
        del info
        gc.collect()
    print "Total #subtomograms:", len(info_copy_json)
    info_json = [_ for _ in info_copy_json if 'pdb_id' in _]
    
    fname = os.path.join(processing_dir_root, opt["file_with_all_subtomograms"])
    print "Writing:", fname
    with open(fname, "w") as f:
        json.dump(info_copy_json, f, indent=2)

    fname = os.path.join(processing_dir_root, opt["file_with_all_subtomograms_that_have_pdb_id"])
    print "Writing:", fname
    with open(fname, "w") as f:
        json.dump(info_json, f, indent=2)

if __name__=='__main__':
    main()
