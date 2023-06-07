#!/usr/bin/env python

# given a subtomogram list file (from stdin) in json format, select the subtomograms with particular cluster labels


from subtomogram_cluster_sel import tomominer.cluster_sel

if __name__ == '__main__':
    import json
    
    import sys
    lbls = [int(sys.argv[i]) for i in range(1, len(sys.argv))]
    lbls = set(lbls)

    conf = json.load(sys.stdin)
    
    conf_f = cluster_sel(conf, lbls, False)
       
    json.dump(conf_f, sys.stdout, indent=2)

