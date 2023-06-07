#!/usr/bin/env python



# functions for packing spheres
import sys
import numpy as N

from multiprocessing import Pool

import tomominer.statistics.util as SU


# given a cubic region and sphere radius and proportions, pack the spheres using dynamics so that the spheres are filled up the region
def pack(op):

    pool = Pool()

    # collect sphere radius and proportions
    rad = {}
    prop = {}
    for r in op['spheres']:
        i = int(r['id'])
        rad[i] = float(r['rad'])
        prop[i] = float(r['proportion']);        assert(prop[i] > 0)


    rad_max = max(rad[_] for _ in rad)
    rad_min = min(rad[_] for _ in rad)

    if False:
        grid_spacing = rad_max
    else:
        grid_spacing = rad_min      # use rad_min so that the initially balls are compact, these balls will expand, the expansion speed will be limited using max_velocity and damping

    #--------------------------------
    # first arrange the particles in to a (expanded) grid with distance to be rad_max (more compact than 2 * rad_max, so that they initially tend to push each other away)

    # generate a grid of ids whos frequency corresponds to the proportion
    container_siz = N.array(op['container']['size'], dtype=N.float)
    grid_siz = N.ceil( container_siz / rad_max);     grid_siz.astype(int)

    ids = generate_id_grid(grid_siz=grid_siz, prop=prop)

    rad_v = N.array( list(rad[_] for _ in ids) )     # vector of radius

    rad_v_enlarge = rad_v


    # generate initial coordinate
    x = coordinate_init(grid_siz=grid_siz, spacing=grid_spacing)
    v = N.zeros(x.shape)        # velocity

    iter_n = 0
    while iter_n < op['iter_max']:
        print '\r', iter_n, '  ',
        x, v = coordinate_evolve(x=x, v=v, rad=rad_v_enlarge, potential_well_depth=rad_max, max_velocity=grid_spacing, pool=pool)
        iter_n += 1

    return {'ids':ids.tolist(), 'x':x.tolist(), 'r':rad_v.tolist()}


    

def generate_id_grid(grid_siz, prop):

    return SU.proportion_sample(size=N.prod(grid_siz), prop=prop)



def coordinate_init(grid_siz, spacing):
    x_grid = N.mgrid[0:grid_siz[0], 0:grid_siz[1], 0:grid_siz[2]]
    x = N.zeros( [N.prod(grid_siz), 3] )
    for dim_i in range(3):      x[:,dim_i] = x_grid[dim_i].flatten()

    x = x * spacing

    return x
 

# evolve particle coordinates acording to force field generated according to Lennard-Jones potential
def coordinate_evolve(x, v, rad, potential_well_depth, max_velocity=None, damping_factor=0.5, time_interval=0.01, pool=None):
    assert( damping_factor <  1.0)       # dampling factor must be less than 1, otherwise positive feedback
    assert(damping_factor >= 0.0)

    damping = damping_factor * (1.0 / time_interval)


    n_chunk = int(   max(1, N.ceil(float(x.shape[0]) / (pool._processes*2)))   )

    inds = list(range(len(rad)))
    pool_res = []
    while len(inds) > 0:

        inds_t = inds[:n_chunk]

        pool_res.append( pool.apply_async(func=calculate_force_field_local, kwds={'x':x, 'v':v, 'rad':rad, 'potential_well_depth':potential_well_depth, 'inds':inds_t})  )

        inds = inds[n_chunk:]


    ff = N.zeros(v.shape)       # force field
    for pool_res_t in pool_res:
        re = pool_res_t.get()
        ff[re['inds'],:] = re['f']

    v_new = v + ff * time_interval
    v_new = v_new - v_new*damping*time_interval     # besides pushing the balls using force field, we also add a damping to make the particles stablize

    
    # thresholding max velocity, this is used to prevent particles to move too fast when force is strong
    v_new_norm = N.sqrt(   N.square(v_new).sum(axis=1)  )
    v_new_norm_ind = (v_new_norm > max_velocity)
    for dim in range(v_new.shape[1]):       v_new[v_new_norm_ind, dim] = (v_new[v_new_norm_ind, dim] / v_new_norm[v_new_norm_ind]) * max_velocity

    x_new = x + (v_new * time_interval)


    ff_mean = N.sqrt( N.square(ff).sum(axis=1) ).mean()
    vd_mean = N.sqrt( N.square(v_new - v).sum(axis=1) ).mean()     # mean velocity change
    v_mean = N.sqrt( N.square(v_new).sum(axis=1) ).mean()     # mean velocity
    
    print 'ff:', ff_mean, '  vd:', vd_mean, '  v', v_mean, '   ', x_new.max(axis=0) - x_new.min(axis=0), ' ', 
    #sys.stdout.flush()

    return (x_new, v_new)
    

def calculate_force_field_local(x, v, rad, potential_well_depth, inds=None):

    f_local = N.zeros( [len(inds), x.shape[1]] )

    for i, ind_i in enumerate(inds):
        x_diff = x - N.tile(x[ind_i,:], [x.shape[0], 1])

        # calculate distance to x_t
        r = N.sqrt( N.square(x_diff).sum(axis=1) )

        r_ind = (r > 0)


        rad_sum = rad + rad[i]
        dp = N.zeros(r.shape)

        if True:
            '''
            calculate differential of Lennard-Jones potential as force (see wiki force field)

            V = \varepsilon \left[\left( \frac {r_{m}} {r} \right)^{12} -  2 \left( \frac {r_{m}} {r} \right)^6 \right]\\
            {d V \over dr} = \varepsilon \left[12\left( \frac {r_{m}} {r} \right)^{11} r_m (-r^{-2}) -  2\times 6 \left( \frac {r_{m}} {r} \right)^5 r_m (-r^{-2}) \right]\\
            = 12\varepsilon \left[\left( \frac {r_{m}} {r} \right)^{11} -  \left( \frac {r_{m}} {r} \right)^5 \right] r_m (-r^{-2}) \\

            V = \varepsilon \left[\left( \frac {r_{m}} {r} \right)^{4} -  2 \left( \frac {r_{m}} {r} \right)^2 \right]\\
            {d V \over dr} = \varepsilon \left[4\left( \frac {r_{m}} {r} \right)^{3} r_m (-r^{-2}) -  2\times 2 \left( \frac {r_{m}} {r} \right) r_m (-r^{-2}) \right]\\
            = 12\varepsilon \left[\left( \frac {r_{m}} {r} \right)^{3} -  \left( \frac {r_{m}} {r} \right) \right] r_m (-r^{-2}) \\
            '''



            dp[r_ind] = 12.0 * potential_well_depth * (N.power(rad_sum[r_ind] / r[r_ind], 11) - N.power(rad_sum[r_ind] / r[r_ind], 5)) * rad_sum[r_ind] / (-N.square(r[r_ind]))

            #dp[r_ind] = dp[r_ind] + 2*r[r_ind]     # append "r^2" term to Lennard-Jones potential so that the particles won't get far away

        else:

            # a simple potential
            dp[r_ind] = - N.square( rad_sum[r_ind] / r[r_ind] ) + 1  # the first terms is in form of Coulomb's law, the second term (1) is a constant attraction force

        # calculate individual forces
        f = N.zeros(x_diff.shape)
        for dim in range(x_diff.shape[1]):
            f[r_ind,dim] = (x_diff[r_ind,dim] / r[r_ind]) * dp[r_ind]

        f_local[i, :] = f.sum(axis=0)
    
    return {'inds':inds, 'f':f_local}








if __name__ == '__main__':
    
    import sys
    op_file = sys.argv[1]
    out_file = sys.argv[2]


    import json
    with open(op_file) as f:       op = json.load(f)

    p = pack(op)

    with open(out_file, 'w') as f:      json.dump(p, f, indent=2)



'''

# example command

~/ln/tomominer/tomominer/geometry/pack/sphere.py ~/ln/tomominer/tomominer/geometry/pack/op-example.json /tmp/out.json


'''




'''

% matlab code for inspection


addpath( strcat( getenv('HOME'), '/ln/frequent_structure/code')  );






clear all;
close all;

d = JSON.loadjson('/media/sf_tmp/out.json');
d = JSON.loadjson('/tmp/out.json');


[x1, x2, x3] = sphere();


r = d.r * 1.0;



hold on
for i = 1 : size(d.x, 1)

    surf(x1*r(i) + d.x(i,1), x2*r(i) + d.x(i,2), x3*r(i) + d.x(i,3));

end

axis equal;
axis tight;






for id_i = unique(d.ids)

    figure;

    hold on
    for i = find(d.ids == id_i)
        

        surf(x1*r(i) + d.x(i,1), x2*r(i) + d.x(i,2), x3*r(i) + d.x(i,3));

    end

    axis equal;
    axis tight;

end






'''





