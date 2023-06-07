#!/usr/bin/env python


# create subfolders for each matched cluster averages, and export subset of subtomograms that matches each cluster averages

if __name__ == '__main__':
    
    import json
    with open('data_config.json') as f:     data = json.load(f)



    import pickle
    with open('cluster_average_select.pickle') as f:    st = pickle.load(f)



    import os
    st = st['selected_templates']
    for c in st:
        data_t = []
        for r in data:
            if str(r['template']['subtomogram']) == str(st[c]['subtomogram']):
                data_t.append(r)

        tem_dir = os.path.join( os.getcwd(), 'tem_%d'%(c) )
        if not os.path.isdir(tem_dir):            os.mkdir(tem_dir)

        with open(os.path.join(tem_dir, 'data_config.json'), 'w') as f:  json.dump(data_t, f, indent=2)

        print c, len(data_t), tem_dir


