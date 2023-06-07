#!/usr/bin/env python


'''
using out distriuted system, parallelly perform template search, written according to FreqStructureFindingByWatershedGrowth.template_search()

~/ln/tomominer/tomominer/template/search/standard_scanning/distributed/run.py

'''



import os
import shutil
import json
import pickle
import copy
import uuid
import numpy as N
import scipy.spatial.distance as SSD

import tomominer.io.file as IOF
import tomominer.geometry.rotate as GR
import tomominer.filter.convolve as FC
import tomominer.filter.normalized_cross_correlation as FNCC
import tomominer.filter.local_extrema as FLE
import tomominer.geometry.point_cloud.redundant_points as GPR
import tomominer.model.util as MU

# parameters: mk: map file name,    tk: template        op: options
def template_match_parallel(op):

    angle_interval = (float(op['angle_interval']) / 180) * N.pi

    # first make a list of all possible rotation angles
    phi_s = N.arange(start=(-N.pi), stop=N.pi, step=angle_interval).tolist()
    theta_s = N.arange(start=0.0, stop=N.pi, step=angle_interval).tolist();         theta_s.append(N.pi)
    psi_s = N.arange(start=(-N.pi), stop=N.pi, step=angle_interval).tolist()

    angles = []
    for phi in phi_s:
        for theta in theta_s:
            for psi in psi_s:
                angles.append(   (phi, theta, psi)   )

    print 'scanning', len(angles), 'angles'

    from tomominer.parallel.queue_master import QueueMaster
    qhost = os.getenv('HOSTNAME')
    qport = 5011
    runner  = QueueMaster(qhost, qport)

    if 'worker_num' not in op:      op['worker_num'] = runner.work_queue.get_worker_number()
    n_chunk=max(op['min_chunk_size'], int( N.ceil(len(angles) / (op['worker_num']+10)) ))

    print 'worker num', op['worker_num'], 'n_chunk', n_chunk

    tasks = []
    angles_t = angles
    while angles_t:
        op_t = {'tmp_dir':op['tmp_dir'], 'template':op['template'], 'map':op['map'], 'mode':op['mode'],'angles':angles_t[:n_chunk], 'task_id':len(tasks)}
        tasks.append(runner.task( module='tomominer.template.search.standard_scanning.run', method='template_match_single_job', kwargs={'op':op_t} ))

        angles_t = angles_t[n_chunk:]

    


    raise Exception('following code need to consider partitions')

    c_max = None
    phi_max = None
    theta_max = None
    psi_max = None

    for r in runner.run__except(tasks):
        r = r.result

        c = N.load(r['c']);        os.remove(r['c'])

        if c_max is None:
            c_max = c

            if 'phi' in r:
                phi_max = N.load(r['phi']);             os.remove(r['phi'])
                theta_max = N.load(r['theta']);             os.remove(r['theta'])
                psi_max = N.load(r['psi']);             os.remove(r['psi'])

        else:

            ind = c > c_max
            c_max[ind] = c[ind]

            phi_max[ind] = N.load(r['phi'])[ind];             os.remove(r['phi'])
            theta_max[ind] = N.load(r['theta'])[ind];             os.remove(r['theta'])
            psi_max[ind] = N.load(r['psi'])[ind];             os.remove(r['psi'])


    print 'finding local maxima'
    x = FLE.local_maxima(c_max)

    c_val = c_max[x]
    phi_val = phi_max[x]
    theta_val = theta_max[x]
    psi_val = psi_max[x]

    x = N.array(x).T

    ps = []
    for i in range(len(c_val)):     ps.append(      {'val':float(c_val[i]), 'x':[_ for _ in x[i,:]], 'angle':[float(phi_val[i]), float(theta_val[i]), float(psi_val[i]) ], 'map_id':op['map_id'], 'template_id':op['template_id']}      )

    # order peaks in ps according to values
    ps = sorted(ps, key=lambda _:(-_['val']))
    for i in range(len(ps)):        ps[i]['id'] = i

    print 'peaks detected', len(ps)

    N.save(os.path.join(op['tmp_dir'], 'match-%d.npy'%(op['match_id'],)), c_max)

    return ps



'''

# given two sets of peaks, if a point in set0 is very close to a point in set1, and if the point
def merge_peaks(p0, p1, r):
    x0 = N.array([_['x'] for _ in p0])
    x1 = N.array([_['x'] for _ in p1])
    d = SSD.cdist(x0, x1)

    f0 = [True] * len(p0)
    f1 = [True] * len(p1)
    for i0 in range(len(p0))
        w = N.where(d[i0,:].flatten() < r)[0].tolist()
        if len(w) == 0:     continue
        for i1 in w:
            if f1[i1] == False: continue

            if p0[i0]['score'] < p1[i1]['score']:
                f0[i0] = False
                break
            else:
                f1[i1] = False
    
    
    p = [p0[_] for _ in range(len(p0)) if f0[_]]
    p.extend( [p1[_] for _ in range(len(p1)) if f1[_]] )

    return p
'''



# when the map is very large, we can partition it, then perform parallel matching for each partition, then collect peaks. Note that we need to filter redundant peaks afterwards
def template_match_parallel__partition(op):

    if not os.path.isdir(op['tmp_dir']):    os.makedirs(op['tmp_dir'])
    
    v = N.load(op['map'])

    import tomominer.image.vol.partition as IVP
    b = IVP.gen_bases(v.shape, nonoverlap_width=op['partition']['nonoverlap_width'], overlap_width=op['partition']['overlap_width'])

    print 'partition size', b.shape
   
    bl = []         # list of bases and other information
    
    for i0 in range(b.shape[0]):
        for i1 in range(b.shape[1]):
            for i2 in range(b.shape[2]):
                bp = N.squeeze(b[i0,i1,i2,:,:])
                partition_dir = os.path.join(op['tmp_dir'], 'partition--%d-%d-%d'%(i0,i1,i2))
                partition_file = os.path.join(partition_dir, 'map.npy')
                bl.append({'dir':partition_dir, 'map_file':partition_file, 'base':bp})
    del b

    with open(os.path.join(op['tmp_dir'], 'partition.pickle'), 'w') as f:       pickle.dump(bl, f)

    # load large map, then store partition files
    vp = None
    for b in bl:
        bp = b['base']
        if not os.path.isdir(b['dir']):            os.makedirs(b['dir'])
        if not os.path.isfile(b['map_file']):
            vp = v[bp[0,0]:bp[0,1], bp[1,0]:bp[1,1], bp[2,0]:bp[2,1]]
            N.save(b['map_file'], vp)

    del vp
    del v

    # search partition one by one and collect peaks
    ps = []
    for b in bl:
        print 'partition dir', b['dir']
        tm_op = copy.deepcopy(op['search'])
        tm_op['map'] = b['map_file']
        tm_op['map_id'] = op['map_id']
        tm_op['tmp_dir'] = b['dir']
        tm_op['match_id'] = op['match_id']

        ps_t = template_match_parallel(tm_op)

        # add base coordinate
        for p in ps_t:
            for i in range(3):
                p['x'][i] += b['base'][i,0]

        ps.extend(ps_t)

    ps = sorted(ps, key=lambda _:(-_['val']))       # order according to score
    rf = GPR.identify(N.array([_['x'] for _ in ps]), op['partition']['redundance_distance'])
    ps = [ps[_] for _ in range(len(ps)) if rf[_]]

    return ps


def template_match_single_job(self, op):

    re = {}
    re['c'] = os.path.join(op['tmp_dir'], '%d-c.npy'%(op['task_id'], ))
    re['phi'] = os.path.join(op['tmp_dir'], '%d-phi.npy'%(op['task_id'], ))
    re['theta'] = os.path.join(op['tmp_dir'], '%d-theta.npy'%(op['task_id'], ))
    re['psi'] = os.path.join(op['tmp_dir'], '%d-psi.npy'%(op['task_id'], ))

    if os.path.isfile(re['c']) and os.path.isfile(re['phi']) and os.path.isfile(re['theta']) and os.path.isfile(re['psi']):     return re


    t = N.load(op['template'])
    tm = N.isfinite(t)      # real space mask
    t_mean = t[tm].mean()
    t[N.logical_not(tm)] = t_mean
    tm = tm.astype(N.float)


    v = N.load(op['map'])

    c_max = None
    phi_max = None
    theta_max = None
    psi_max = None

    for (phi, theta, psi) in op['angles']:
        tr = GR.rotate(t, angle=(phi, theta, psi), default_val=t_mean)

        if op['mode'] == 'convolve':
            c = FC.convolve(v=v, t=tr)
        elif op['mode'] == 'normalized-cor':
            tmr = GR.rotate(tm, angle=(phi, theta, psi), default_val=0.0)
            tr[tmr < 0.5] = float('NaN')
            c = FNCC.cor(v=v, t=tr)
        else:
            raise Exception('mode')

        if c_max is None:
            c_max = c
            phi_max = N.zeros(c.shape) + phi
            theta_max = N.zeros(c.shape) + theta
            psi_max = N.zeros(c.shape) + psi

        else:

            ind = (c > c_max)
            c_max[ind] = c[ind]
            phi_max[ind] = phi
            theta_max[ind] = theta
            psi_max[ind] = psi

   
    if not os.path.isfile(re['c']):    
        N.save(re['c'], c_max)
        N.save(re['phi'], phi_max)
        N.save(re['theta'], theta_max)
        N.save(re['psi'], psi_max)

    return re

    


    


def main():    
    with open('template_search__op.json') as f:     op = json.load(f)

    op['tmp_dir'] = os.path.join(op['tmp_dir'], str(uuid.uuid4()))
    print 'tmp dir', op['tmp_dir']
    
    if os.path.isdir(op['tmp_dir']):        shutil.rmtree(op['tmp_dir'])
    os.makedirs(op['tmp_dir'])

    re = []
    for map_t in op['maps']:
        # first read whole map, and store it into the temporary folder in numpy array format
        v = IOF.read_mrc(map_t['vol_file'])['value']

        tmp_map = os.path.join(op['tmp_dir'], 'map-%d.npy'%(map_t['id'],))
        N.save(tmp_map, v)

        for tm in op['templates']:
            t = IOF.read_mrc(tm['vol_file'])['value']
            if op['inverse_template']:  t = -t

            # apply a spherical mask. In future, can use segmentation to obtain a more precise mask
            t[MU.sphere_mask(t.shape) == 0] = float('NaN')

            tmp_template = os.path.join(op['tmp_dir'], 'template-%d.npy'%(tm['id'],))
            N.save(tmp_template, t)

            if 'partition' in op:
                tmp_op = {}
                tmp_op['partition'] = copy.deepcopy(op['partition'])
                tmp_op['search'] = copy.deepcopy(op['search'])
                tmp_op['search']['template'] = tmp_template
                tmp_op['search']['template_id'] = tm['id']
                tmp_op['map'] = tmp_map
                tmp_op['map_id'] = map_t['id']
                tmp_op['match_id'] = len(re)
                tmp_op['tmp_dir'] = os.path.join(op['tmp_dir'], 'match-%d'%(tmp_op['match_id'],))

                tm_re = template_match_parallel__partition(tmp_op)
            else:
                tm_op = copy.deepcopy(op['search'])
                tm_op['template'] = tmp_template
                tm_op['template_id'] = tm['id']
                tm_op['map'] = tmp_map
                tm_op['map_id'] = map_t['id']
                tm_op['tmp_dir'] = op['tmp_dir']
                tm_op['match_id'] = len(re)

                tm_re = template_match_parallel(tm_op)

            re.extend(tm_re)


    with open('template_search__out.json', 'w') as f:   json.dump(re, f, indent=2)


if __name__ == "__main__":
    main()




'''

# test code, first generate a volume

import tomominer.model.util as MU
t = MU.generate_toy_model(dim_size=32)

import numpy as N
import tomominer.geometry.rotate as GR

ealarge_ratio = 4
v0 = GR.rotate(t, c2=N.array(t.shape) * 1, siz2=N.array(t.shape)*ealarge_ratio, default_val=0.0)
v1 = GR.rotate(t, angle=N.array([N.pi/2,0,0]), c2=N.array(t.shape) * 3, siz2=N.array(t.shape)*ealarge_ratio, default_val=0.0)
v = v0+v1

import matplotlib
matplotlib.use('Qt4Agg')

import tomominer.image.vol.util as IVU
IVU.dsp_cub(v)


import tomominer.io.file as IF
IF.put_mrc(t, '/tmp/template.mrc', overwrite=True)
IF.put_mrc(v, '/tmp/map.mrc', overwrite=True)


'''


'''
template_search__op.json

{

    "program" : "~/ln/tomominer/tomominer/template/search/standard_scanning/run.py",

    "maps" : [
        {
            "id" : 0,
            "vol_file" : "/tmp/map.mrc"
        }
    ],

    "templates" : [
        {
            "id" : 0,
            "vol_file" : "/tmp/template.mrc"
        }
    ],

    "tmp_dir" : "/home/rcf-47/mxu/tmp-staging/template-search",

    "inverse_template" : false,

    "search" : {
        "angle_interval" : 90,
        "min_chunk_size" : 2,
        "worker_num": 10,
        "mode" : "convolve",
        "mode bak" : "normalized-cor"
    },

    "partition" : {
        "nonoverlap_width" : 64,
        "overlap_width" : 32,
        "redundance_distance" : 5.0
    }

}




'''


'''

# plot c_max, maximum cross correlation map

import tomominer.io.file as IF
v = IF.read_mrc('/tmp/map.mrc')['value']

import numpy as N
c = N.load('/home/rcf-47/mxu/tmp-staging/template-search/match-0/partition--0-0-0/match-0.npy')


import matplotlib
matplotlib.use('Qt4Agg')

import tomominer.image.vol.util as IVU
IVU.dsp_cub(v)
IVU.dsp_cub(c)


'''

