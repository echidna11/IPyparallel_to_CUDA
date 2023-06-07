#!/usr/bin/env python


# calculate contingency table of between cluster labels and template matches



def stat(dj):

    from collections import defaultdict
    c = defaultdict(dict)

    for d in dj:
        if 'cluster_label' not in d:     continue
        if 'template' not in d:     continue

        l = int(d['cluster_label'])

        t = str(d['template']['subtomogram'])

        if l not in c[t]:   c[t][l] = 0
        c[t][l] += 1

    return c


def stat_print(f, c):
   
    # prepare and write header
     
    labels = set()
    for t in c:
        for l in c[t]:
            labels.add(l)

    labels = list(labels)
    labels.sort()

    for l in labels:    
        f.write('\t')
        f.write(repr(l))
    f.write('\n')

    # write table body
    for t in c:
        f.write(t)

        for l in labels:
            f.write('\t')
            if l in c[t]:
                f.write(repr(c[t][l]))
            else:
                f.write(repr(0))
        f.write('\n')




if __name__ == '__main__':
    

    import os
    import json

    # walk through every subdir, find and record last passes
    for root, sub_folders, files in os.walk(os.getcwd()):

        data_config_file = os.path.join(root, 'data_config.json')

        if not os.path.exists(data_config_file):        continue

        with open(data_config_file) as f:   dj = json.load(f)

        c = stat(dj)

        if len(c) == 0:    continue

        print root

        with open(os.path.join(root, 'cluster_consistency_with_match.txt'), 'w')  as f:       stat_print(f, c)


