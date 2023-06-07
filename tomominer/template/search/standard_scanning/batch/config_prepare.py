#!/usr/bin/env python



'''
prepare configuration on rotation angle partition and tomogram partition


~/ln/tomominer/tomominer/template/search/standard_scanning/batch/config_prepare.py

'''


import os, json, uuid
import cPickle as pickle
import numpy as N

import tomominer.io.file as IF

def config_prepare(op):

    op['out_dir'] = os.path.abspath(op['out_dir'])
    if not os.path.isdir(op['out_dir']):    os.makedirs(op['out_dir'])

    op['this_file'] = os.path.abspath(op['this_file'])

    mrc = IF.read_mrc_header(op['tomogram'])['MRC']

    bl = partition_list(shape=[mrc['nx'], mrc['ny'], mrc['nz']], nonoverlap_width=op['partition']['nonoverlap_width'], overlap_width=op['partition']['overlap_width'], out_dir=os.path.join(op['out_dir'], 'partitions'))

    angles = angle_list( (float(op['angle_interval']) / 180) * N.pi )

    jobs = []
    
    for b in bl:
        angles_t = angles
        while angles_t:
            jobs.append({'angles':angles_t[:op['angle_chunk_size']], 'partition':b, 'id':str(uuid.uuid4())})
            angles_t = angles_t[op['angle_chunk_size']:]

    print len(jobs), 'jobs generated'

    
    jobs_file = os.path.join(op['out_dir'], 'jobs.pickle')
    with open(jobs_file, 'wb') as f:      pickle.dump({'jobs':jobs, 'angles':angles, 'partitions':bl}, f, protocol=-1)

    job_file_list = []
    for jobs_t in jobs:
        tid = jobs_t['id']
        job_dir = os.path.join(op['out_dir'], 'jobs', tid)
        if not os.path.isdir(job_dir):     os.makedirs(job_dir)

        job_file = os.path.join(job_dir, 'job.pickle')
        with open(job_file, 'wb') as f:        pickle.dump({'id':tid, 'job':jobs_t, 'job_dir':job_dir, 'out_dir':os.path.join(job_dir, 'out'), 'stat_out':os.path.join(job_dir, 'stat.json'), 'this_file':job_file, 'config_file':op['this_file'], 'config_stat_out_file':op['config_stat_out_file']}, f, protocol=-1)
        
        job_file_list.append({'id':tid, 'job_file':job_file})

    return {'jobs_file':jobs_file, 'job_file_list':job_file_list, 'template_tmp':os.path.join(op['out_dir'], 'template.npy')}




def angle_list(angle_interval):

    # first make a list of all possible rotation angles
    phi_s = N.arange(start=(-N.pi), stop=N.pi, step=angle_interval).tolist()
    theta_s = N.arange(start=0.0, stop=N.pi, step=angle_interval).tolist();         theta_s.append(N.pi)
    psi_s = N.arange(start=(-N.pi), stop=N.pi, step=angle_interval).tolist()

    angles = []
    for phi in phi_s:
        for theta in theta_s:
            for psi in psi_s:
                angles.append(   (phi, theta, psi)   )

    return angles


# when the map is very large, we can partition it, then perform parallel matching for each partition, then collect peaks. Note that we need to filter redundant peaks afterwards
def partition_list(shape, nonoverlap_width, overlap_width, out_dir):

    import tomominer.image.vol.partition as IVP
    b = IVP.gen_bases(shape, nonoverlap_width=nonoverlap_width, overlap_width=overlap_width)

    print 'partition size', b.shape
   
    bl = []         # list of bases and other information
    
    for i0 in range(b.shape[0]):
        for i1 in range(b.shape[1]):
            for i2 in range(b.shape[2]):
                bp = N.squeeze(b[i0,i1,i2,:,:])
                partition_dir = os.path.join(out_dir, 'partition--%d-%d-%d'%(i0,i1,i2))
                partition_file = os.path.join(partition_dir, 'map.npy')
                bl.append({'dir':partition_dir, 'map_file':partition_file, 'base':bp, 'id':str(uuid.uuid4())})

    return bl



def main():
    config_file = os.path.abspath('config_prepare__op.json')
    with open(config_file) as f:     op = json.load(f)
    op['this_file'] = config_file
    op['config_stat_out_file'] = os.path.abspath(op['config_stat_out_file'])

    c = config_prepare(op)

    with open(op['config_stat_out_file'], 'w') as f:        json.dump(c, f, indent=2)


if __name__ == '__main__':
    main()



