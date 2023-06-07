

def cluster_sel(conf, lbls, include=True):
    conf_f = []
    for record in conf:
        if 'cluster_label' not in record: print record; raise Exception
        
        if include:   # in this case, want to include selected clusters 
            if not int(record['cluster_label']) in lbls: continue
        else:       # in this case, want to exclude selected clusters
            if int(record['cluster_label']) in lbls: continue
                        
        conf_f.append(record)
 
    return conf_f

