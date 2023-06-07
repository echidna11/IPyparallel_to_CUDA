#!/usr/bin/env python


'''
given a collection of density maps, generate templates


~/ln/tomominer/tomominer/template/library/generate_templates.py

'''


import os, json, copy

import tomominer.template.generate.map_to_template as TGM

def main():
    with open('generate_templates__op.json') as f:      op = json.load(f)
    op['out dir'] = os.path.abspath(op['out dir'])

    if not os.path.isdir(op['out dir']):    os.makedirs(op['out dir'])

    with open(op['map op']) as f:       mop = json.load(f)
    with open(op['map in']) as f:       ps = json.load(f)


    fs = []
    for pid in ps:
        cop = copy.deepcopy(op['convert'])
        cop['ctf']['pix_size'] = mop['situs']['spacing']
        cop['map_file'] = ps[pid]
        cop['template_file_out'] = os.path.join(op['out dir'], pid + '-vol.mrc')
        cop['mask_file_out'] = os.path.join(op['out dir'], pid + '-mask.mrc')
        TGM.convert(cop)

        fs.append({'pid':pid, 'subtomogram':cop['template_file_out'], 'mask':cop['mask_file_out']})

    with open(op['stat out'], 'w') as f:     json.dump(fs, f, indent=2)

if __name__ == '__main__':
    main()



'''
related code
~/ln/tomominer/tomominer/template/generate/map_to_template.py


the output can be used by 
~/ln/tomominer/tomominer/template/search/fast_align/reference/multi.py

'''

