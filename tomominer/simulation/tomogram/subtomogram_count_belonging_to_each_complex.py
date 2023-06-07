#code to count no. of instances actually in picked particles

import json, os
def main():
    with open("model_tomogram_subtomogram__config.json") as f:    op = json.load(f)
    with open(op["model_generation_config_file"]) as f:    model_op = json.load(f)
    indir = os.path.join(op['out_dir'], str(model_op['copy_number']['total']), '%d-%d-%d'%(model_op['packing']['param']['box']['x'], model_op['packing']['param']['box']['y'], model_op['packing']['param']['box']['z']))
    with open(os.path.join(indir, "info.json")) as f:    info = json.load(f)

    infe = [_ for _ in info if 'pdb_id' in _]
    pdbs = set([_['pdb_id'] for _ in infe])

    print "Number of instances in info.json:", len(info)
    print "Number of instances with identified pdb_id:", len(infe)
    print "Number of PDBs identified:", len(pdbs)

    count = {}
    tomogram_ids = list(set([_['tomogram_id'] for _ in infe]))
    for t_id in tomogram_ids:
	count[t_id] = {}
    	for i in pdbs:
        	count[t_id][i] = sum([1 for _ in infe if _['pdb_id']==str(i) and _['tomogram_id']==t_id])

    out_file = "subtomogram_count_belonging_to_each_complex.json"
    with open(os.path.join("./", out_file), 'w') as f:  json.dump(count, f, indent = 2)

if __name__ == '__main__':
    main()
