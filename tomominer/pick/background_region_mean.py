#!/usr/bin/env python



# given a manually defined empty background region, calculate its mean


if __name__ == '__main__':

    import json
    with open('background_region_mean__config.json') as f:      op = json.load(f)

    with open(op['tomogram_info_file']) as f:     tom_info = json.load(f)
    tom_info = {_['id']:_ for _ in tom_info}

    import tomominer.io.file as IF
    import numpy as N

    re = []
    for t in op['tomograms']:
        print 'tomogram', t['tomogram_id']
        v = IF.read_mrc(tom_info[t['tomogram_id']]['vol_file'], show_progress=True)['value']

        # clip the map according to given corners
        c = N.array(t['corner'])
        v = v[c[0,0]:c[0,1], c[1,0]:c[1,1], c[2,0]:c[2,1]]

        re.append(  {'tomogram_id':t['tomogram_id'], 'mean':float(v.mean())}  )


    with open('background_region_mean__out.json', 'w') as f:     json.dump(re, f, indent=2)



