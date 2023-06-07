import copy
import math
import numpy as N
from numpy.fft import fftn, ifftn, fftshift, ifftshift
import scipy.linalg as SL


import tomominer.image.vol.util as IMU
import tomominer.geometry.ang_loc as GA
import tomominer.geometry.rotate as GR
import tomominer.statistics.regression.linear as SRL

import tomominer.io.file as IF

import tomominer.core as core


def foster_transform_mask(m1, m2, ang, mask_threshold=0.1):

    mid_co = IMU.fft_mid_co(m2.shape)
    if True:
        m2r = GR.rotate(m2, c1=mid_co, c2=mid_co, angle=ang, default_val=0.0)
    else:
        m2r = core.rotate_mask(m2, ang)

    m2r = N.abs(m2r)        # sometimes rotation interpolation brings negative values...

    
    m = m1 * m2r
    mid_co = mid_co.astype(N.int)
    m[mid_co[0], mid_co[1], mid_co[2]] = 1.0        # make sure that the zero frequency component is kept

    return m


def foster_transform_single_vol(v, m):

    vf = fftshift(fftn(v))
    vfi = N.real( ifftshift(ifftn(vf * m)) )

    vfi_f = vfi - N.mean(vfi)
    vfi_f = vfi_f / math.sqrt( N.square(vfi_f).sum() )

    return vfi_f



# gradient based refinement, single step
# see SphericalHarmonicsUtil.gradient_search_map_increment()
def refine(op):
    v1 = op['v1']
    m1 = op['m1']
    v2 = op['v2']
    m2 = op['m2']
    ang = op['ang']
    loc_r = op['loc_r']

    if 'lambda' in op:
        lambda_t = op['lambda']
    else:
        lambda_t = 1.0
    
    ang_epsilon = op['ang_epsilon']
    loc_epsilon = op['loc_epsilon']


    (msk, v1f, v2r, v2rf, dist_sum, cor_t) = refine__calculate_transform(v1=v1, m1=m1, v2=v2, m2=m2, ang=ang, loc_r=loc_r)

    eps = N.finfo(N.float64).eps      # machine epsilon

    J = N.zeros([v1.size, 6])

    for dim_i in range(3):

        if abs(ang_epsilon[dim_i]) > eps:
            ang_epslon_t = N.zeros(3)
            ang_epslon_t[dim_i] = ang_epsilon[dim_i]

            v2rt_p = GR.rotate_pad_mean(v2, angle=(ang + ang_epslon_t), loc_r=loc_r)
            v2rt_m = GR.rotate_pad_mean(v2, angle=(ang - ang_epslon_t), loc_r=loc_r)

            v2rt_pf = foster_transform_single_vol(v2rt_p, msk)
            v2rt_mf = foster_transform_single_vol(v2rt_m, msk)

            J[:, dim_i] = (v2rt_pf.flatten() - v2rt_mf.flatten()) / (2*ang_epsilon[dim_i])              # two sided gradient computation is numberically better. See 12:39 of http://helper.ipam.ucla.edu/wowzavideo.aspx?vfn=10741.mp4&vfd=gss2012

            del v2rt_p, v2rt_m, v2rt_pf, v2rt_mf

        else:
            J[:, dim_i] = 0.0
            
    for dim_i in range(3):
        if abs(loc_epsilon[dim_i]) > eps:
            loc_epslon_t = N.zeros(3)
            loc_epslon_t[dim_i] = loc_epsilon[dim_i]

            v2rt_p = GR.rotate_pad_mean(v2, angle=ang, loc_r=(loc_r + loc_epslon_t))
            v2rt_m = GR.rotate_pad_mean(v2, angle=ang, loc_r=(loc_r - loc_epslon_t))


            v2rt_pf = foster_transform_single_vol(v2rt_p, msk)
            v2rt_mf = foster_transform_single_vol(v2rt_m, msk)


            J[:, dim_i + 3] = (v2rt_pf.flatten() - v2rt_mf.flatten()) / (2*loc_epsilon[dim_i])

            del v2rt_p, v2rt_m, v2rt_pf, v2rt_mf

        else:
            J[:, dim_i + 3] = 0.0
            
                
    JQ = J.transpose().dot(J)
    JQ_t = JQ + lambda_t * N.diag(N.diag(JQ))

    JQ_col_sum = (abs(JQ_t) + abs(JQ_t.transpose())).sum(axis=1)        # for the sum, axis=0 sums along first axis, axis=1 sums along second axis
    non_singular_ind = JQ_col_sum > abs(JQ).max() * 1e-10

    if non_singular_ind.sum() == 0:
        print 'warning: all columns are singular   '
        return None

    if not N.all(non_singular_ind):
        print 'warning: some columns are singular   '


    '''
    if numpy.linalg.cond(JQ_t[non_singular_ind, non_singular_ind]) < 1e-7:
        print JQ_col_sum[non_singular_ind]
        print JQ[non_singular_ind, non_singular_ind]
    '''

    f_dif = v1f.flatten() - v2rf.flatten()


    solve_A = JQ_t[non_singular_ind, :][:, non_singular_ind]        # this extraction of submatrix is in a different syntax to matlab
    solve_b = (J[:, non_singular_ind].transpose()).dot(f_dif)
    #print non_singular_ind
    #print solve_A
    #print solve_b

    param_d = N.zeros(6)
    if True:
        param_d[non_singular_ind] = SL.solve( solve_A,  solve_b, sym_pos=True)
    else:
        param_d[non_singular_ind], _, _, _ = N.linalg.lstsq( solve_A,  solve_b )        # this option does not work well
    #print solve_A, solve_b
    #print param_d

    ang_d = param_d[0:3]
    loc_r_d = param_d[3:6]

    #print ang_d
    #print loc_r_d

    return {'ang_d':ang_d, 'loc_r_d':loc_r_d, 'dist_sum':dist_sum, 'cor':cor_t}


def refine__calculate_transform(v1, m1, v2, m2, ang, loc_r):
    
    msk = foster_transform_mask(m1, m2, ang)
    v1f = foster_transform_single_vol(v1, msk)

    if True:
        v2r = GR.rotate_pad_mean(v2, angle=ang, loc_r=loc_r)
    else:
        v2r = core.rotate_vol_pad_mean(v2, ang, loc_r)

    v2rf = foster_transform_single_vol(v2r, msk)


    # calculate squared distancecor_t
    dist_sum = math.sqrt(N.square(v1f - v2rf).sum())

    # calculate cross correlation
    cor_t = (v1f * v2rf).sum()

    if cor_t > 0.99:    raise Exception('score too good to be true, possible checkboard after transform!!!')

    return (msk, v1f, v2r, v2rf, dist_sum, cor_t)



def gradient_search(op):
    
    ang_epsilon = N.ones(3) * (math.pi * (1.0/180))
    loc_epsilon = N.ones(3) * 1.0

    ang = op['ang']
    loc_r = op['loc_r']

    init_dist_sum = None
    old_dist_sum = None

    for i in range(op['max_iter_num']):
        re = refine({'v1':op['v1'], 'm1':op['m1'], 'v2':op['v2'], 'm2':op['m2'], 'ang':ang, 'loc_r':loc_r, 'ang_epsilon':ang_epsilon, 'loc_epsilon':loc_epsilon})
        ang += re['ang_d']
        loc_r += re['loc_r_d']
        print re['dist_sum'], '\t', (2-re['dist_sum']**2) / 2.0, re['ang_d'], '\t', re['loc_r_d'], '\t', ((GA.angle_zyz_translation_difference(ang1=op['truth']['ang'], loc1_r=op['truth']['loc_r'], ang2=ang, loc2_r=loc_r)) if 'truth' in op else ''), '\t', (GA.angle_difference_zyz(ang1=op['truth']['ang'], ang2=ang) if 'truth' in op else ''), '\t', (N.linalg.norm(loc_r - op['truth']['loc_r']) if 'truth' in op else '')


        #print '--------------------------------------------------'

        if True:
            #raise Exception('This is not a good idea, it will lead to stuck at not close to optimal. Is this true???')

            # adjust epsilons for numerical stability
            ang_epsilon_t = abs(re['ang_d']) / 1
            ind = (ang_epsilon_t < ang_epsilon)
            if sum(ind) > 0:            ang_epsilon[ind] = ang_epsilon_t[ind]

            loc_epsilon_t = abs(re['loc_r_d']) / 1
            ind = (loc_epsilon_t < loc_epsilon)
            if sum(ind) > 0:            loc_epsilon[ind] = loc_epsilon_t[ind]

        if init_dist_sum is None:
            init_dist_sum = re['dist_sum']
            continue

        if old_dist_sum is None:
            old_dist_sum = re['dist_sum']
            continue
            
        if abs(re['dist_sum'] - old_dist_sum) < ( abs(init_dist_sum) *  op['stopping_tolerance']):            break
        
        old_dist_sum = re['dist_sum']          
                
        
    return {'ang':ang, 'loc_r':loc_r, 'dist_sum':re['dist_sum']}                


'''
# testing function, load one toy tomogram, then randomly assign an initial rotation and translation, see if iteratively applying refine() can converge

v1f = '/panfs/cmb-panasas2/mxu/proj/imaging/electron/frequent_structure/data/out/method/classification-system/pytom/thermosomes/particle/test_SNR1wedge30/test_SNR1wedge30_100.mrc'
m1f = '/panfs/cmb-panasas2/mxu/proj/imaging/electron/frequent_structure/data/out/method/classification-system/pytom/thermosomes/particle/test_SNR1wedge30/wedge-mask.mrc'
v2f = '/panfs/cmb-panasas2/mxu/proj/imaging/electron/frequent_structure/data/out/method/classification-system/pytom/thermosomes/particle/test_SNR1wedge30/test_SNR1wedge30_100.mrc'
m2f = '/panfs/cmb-panasas2/mxu/proj/imaging/electron/frequent_structure/data/out/method/classification-system/pytom/thermosomes/particle/test_SNR1wedge30/wedge-mask.mrc'

import numpy as N
import tomominer.io.file as IF
v1 = IF.read_mrc_vol(v1f).astype(N.float)
m1 = IF.read_mrc_vol(m1f).astype(N.float)

ang = N.random.rand(3) * N.pi * 2
loc_r = N.zeros(3)

if False:
    v2 = IF.read_mrc_vol(v2f).astype(N.float)
    m2 = IF.read_mrc_vol(m2f).astype(N.float)
else:
    v2 = N.copy(v1)
    m2 = N.copy(m1)
    ang *= 0.01


import tomominer.align.refine.gradient_refine as ARG
re = ARG.gradient_search({'v1':v1, 'm1':m1, 'v2':v2, 'm2':m2, 'ang':ang, 'loc_r':loc_r, 'max_iter_num':100, 'stopping_tolerance':1e-10})

'''




'''
# testing function, load one toy tomogram, then randomly assign an initial rotation and translation, then make subtomograms, see if iteratively applying refine() can converge


import tomominer.align.fast.util as AFU
tp = AFU.generate_toy_subtomogram_pair()

import numpy as N
import tomominer.align.refine.gradient_refine as ARG
re = ARG.gradient_search({'v1':tp['v1b'], 'm1':tp['m1'], 'v2':tp['v2b'], 'm2':tp['m2'], 'ang':tp['ang'] + N.ones(3)*N.pi*0.05, 'loc_r':tp['loc'] + (N.random.random(3)-0.5)*6.0, 'truth':{'ang':tp['ang'], 'loc_r':tp['loc']}, 'max_iter_num':100, 'stopping_tolerance':1e-10, 'verbose':True})

print re

import tomominer.geometry.ang_loc as GA
print GA.angle_zyz_translation_difference(ang1=tp['ang'], loc1_r=tp['loc'], ang2=re['ang'], loc2_r=re['loc_r'])
print GA.angle_difference_zyz(ang1=tp['ang'], ang2=re['ang'])
print re['loc_r'] - tp['loc']


'''



# stocastic parallel refinement
# see SphericalHarmonicsUtil.gradient_search_parallel()

def stocastic_parallel_search(op):
   
    if 'verbose' not in op:     op['verbose'] = False
    
    angs_org = op['angs_org']
    locs_r_org = op['locs_r_org']
    v1 = op['v1']
    m1 = op['m1'];      assert(m1.shape == v1.shape)
    v2 = op['v2'];    assert(v2.shape == v1.shape)
    m2 = op['m2'];    assert(m2.shape == v1.shape)

    min_dist_record_size = op['min_dist_record_size']
    stopping_tolerance = op['stopping_tolerance']
    lambda_t = op['lambda']



    candidate_num = len(angs_org)

    angs = angs_org
    locs_r = locs_r_org


    
    epsilon_factor = 1

    angs_epsilon = N.ones([candidate_num, 3]) * (N.pi * (1.0/180))
    locs_r_epsilon = N.ones([candidate_num, 3])

    if not op.has_key('dists_sum'):
        # initialize distance scores
        initial_dists_sum = N.zeros(candidate_num)
        for ang_i in range(candidate_num):
            (msk__, v1f__, v2r__, v2rf__, initial_dists_sum[ang_i], cor_t__) = refine__calculate_transform(v1, m1, v2, m2, angs[ang_i], locs_r[ang_i])
    else:
        #print 'dists_sum given'
        initial_dists_sum = op['dists_sum']        



    # random sample according to distance scores, and perform refinement
    dists_sum = N.copy(initial_dists_sum)
    initial_dists_sum_min = max(min(initial_dists_sum), 1e-5)
    old_dists_sum = N.zeros(candidate_num) + float('inf')

    dist_sum_s = N.array([])
    min_dist_sum_s = N.array([])
    min_dists_sum_slope = float('inf')

    if op.has_key('top_num'):
        # we consider only a number of top hits for refinement
        top_num = min(op['top_num'], candidate_num)
    else:
        top_num = candidate_num



    initial_dists_sum_sort = N.sort(initial_dists_sum);        initial_dists_sum_sort_cut = max(initial_dists_sum_sort[:top_num])
    search_flags = N.zeros(candidate_num);      search_flags[initial_dists_sum <= initial_dists_sum_sort_cut] = True
            
    pose_difference_threshold = 1.0         # this number is very useful for decreasing the number of search candidates, but may ignore some better candidates, use wisely
        
    if not op.has_key('max_iteration_num'):
        op['max_iteration_num'] = 1000
    
    iteration_num = 0
    
    while iteration_num < op['max_iteration_num']:
       
        # probability sampling
        inds = N.flatnonzero(search_flags == True)
        minimum_scale = math.pow(10.0, (10.0 / float(len(inds) - 1)))       # important! need to use float representation, otherwise integer division lead such minimum_scale mechanism useless !!!
        dists_sum_i = (-dists_sum[inds]).argsort()      # use minus to find indics of descenting order
        inds = inds[dists_sum_i];

        if True:
            ## for code verification
            for ind_i in range(len(inds)-1):
                assert(search_flags[inds[ind_i]] == True)                       # make sure it is a candidate
                assert(dists_sum[inds[ind_i]] >= dists_sum[inds[ind_i+1]])      # make sure descending order


        probs_t = N.zeros(len(inds));      probs_t[0] = 1.0;
        for i in range(1, len(inds)):
            probs_t[i] = probs_t[i-1] * max(minimum_scale, dists_sum[inds[i-1]] / N.max((dists_sum[inds[i]], 1e-5)))

        probs = N.zeros(candidate_num);             probs[inds] = probs_t / sum(probs_t)
       
        # calculate accumulative sum
        probs_sum = N.zeros(len(probs));     prob_sum = 0.0;
        for prob_i in range(len(probs)):
            prob_sum += probs[prob_i]
            probs_sum[prob_i] = prob_sum
        del prob_sum

        while True:
            ang_i = min(N.flatnonzero(probs_sum > N.random.rand()))
            if search_flags[ang_i]:     break

        # perform a single refinement step of the selected rigid transform candidate, denoted by ang_i
        re = refine({'v1':v1, 'm1':m1, 'v2':v2, 'm2':m2, 'ang':angs[ang_i], 'loc_r':locs_r[ang_i], 'ang_epsilon':angs_epsilon[ang_i], 'loc_epsilon':locs_r_epsilon[ang_i], 'lambda':lambda_t})

        if re is None:            break

        dists_sum[ang_i] = re['dist_sum']      ;       assert  N.isfinite(re['dist_sum'])
        angs[ang_i] += re['ang_d']
        locs_r[ang_i] += re['loc_r_d']

        if True:
            #raise Exception('This is not a good idea, it will lead to stuck at not close to optimal. is this true???')

            # adjust epsilons for numerical stability
            ang_epsilon_t = abs(re['ang_d']) / epsilon_factor 
            ind = (ang_epsilon_t < angs_epsilon[ang_i])
            if sum(ind) > 0:            angs_epsilon[ang_i][ind] = ang_epsilon_t[ind]

            loc_epsilon_t = abs(re['loc_r_d']) / epsilon_factor
            ind = (loc_epsilon_t < locs_r_epsilon[ang_i])
            if sum(ind) > 0:            locs_r_epsilon[ang_i][ind] = loc_epsilon_t[ind]

        
        # if two transforms are very close, delete one of them
        ang_i1 = ang_i
        for ang_i2 in range(candidate_num):
            if ang_i2 == ang_i1:
                continue
            
            if search_flags[ang_i2] == False:
                continue

            dif = GA.angle_zyz_translation_difference(angs[ang_i1], locs_r[ang_i1], angs[ang_i2], locs_r[ang_i2]);

            if dif < pose_difference_threshold:
                if dists_sum[ang_i1] < dists_sum[ang_i2]:
                    search_flags[ang_i2] = False
                    continue
                else:
                    search_flags[ang_i1] = False;
                    break

        if op['verbose']:        
            print 'ang_i ' + repr(ang_i),       
            print '    search_flags ' + repr(int(sum(search_flags == True))),
        
        # to decide stop
        dists_sum_current = dists_sum[ang_i];       dist_sum_s = N.append(dist_sum_s, dists_sum_current)
        dists_sum_min = N.min(dists_sum);         min_dist_sum_s = N.append(min_dist_sum_s, dists_sum_min)

        if op['verbose']:
            print '    dists_sum_min ' + repr(dists_sum_min),
            print '    dists_sum_current ' + repr(dists_sum_current),
        
        if (iteration_num > min_dist_record_size) and (dists_sum_current == dists_sum_min):
            iteration_ind_t = N.flatnonzero(abs(dist_sum_s - min_dist_sum_s) == 0);       iteration_ind_t = abs( N.sort((-iteration_ind_t)) );    # use minus for descent order sort

            if len(iteration_ind_t) >= min_dist_record_size:
                min_dist_sum_s_t = min_dist_sum_s[iteration_ind_t[:min_dist_record_size]]
                min_dist_sum_s_t_rev = list(min_dist_sum_s_t);      min_dist_sum_s_t_rev.reverse()          # reverse order for regression
                min_dist_sum_s_t_rev = N.array(min_dist_sum_s_t_rev)

                min_dists_sum_slope = SRL.one_var(N.array(range(len(min_dist_sum_s_t))), (min_dist_sum_s_t_rev / initial_dists_sum_min))[0]

                if op['verbose']:      print '      min_dists_sum_slope ' + repr(min_dists_sum_slope),

                if abs(min_dists_sum_slope) < stopping_tolerance:                    break

        if op['verbose'] and ('truth' in op):        print '   diff to truth', (GA.angle_zyz_translation_difference(ang1=op['truth']['ang'], loc1_r=op['truth']['loc_r'], ang2=angs[ang_i], loc2_r=locs_r[ang_i])), '\t', GA.angle_difference_zyz(ang1=op['truth']['ang'], ang2=angs[ang_i]), '\t', N.linalg.norm(locs_r[ang_i] - op['truth']['loc_r']),

    
        old_dists_sum[ang_i] = dists_sum[ang_i]
        
        iteration_num = iteration_num + 1

        if op['verbose'] : print

    cors = (2 - N.square(dists_sum)) / 2.0

    dists_sum_i = dists_sum.argsort()
    dists_sum_i = dists_sum_i[0];           assert(dists_sum[dists_sum_i] == min(dists_sum))

    return {'dist_sum':dists_sum[dists_sum_i], 'cor':cors[dists_sum_i], 'ang':angs[dists_sum_i], 'loc_r':locs_r[dists_sum_i]}



'''
# test code

v1f = '/panfs/cmb-panasas2/mxu/proj/imaging/electron/frequent_structure/data/out/method/classification-system/pytom/thermosomes/particle/test_SNR1wedge30/test_SNR1wedge30_100.mrc'
m1f = '/panfs/cmb-panasas2/mxu/proj/imaging/electron/frequent_structure/data/out/method/classification-system/pytom/thermosomes/particle/test_SNR1wedge30/wedge-mask.mrc'
v2f = '/panfs/cmb-panasas2/mxu/proj/imaging/electron/frequent_structure/data/out/method/classification-system/pytom/thermosomes/particle/test_SNR1wedge30/test_SNR1wedge30_101.mrc'
m2f = '/panfs/cmb-panasas2/mxu/proj/imaging/electron/frequent_structure/data/out/method/classification-system/pytom/thermosomes/particle/test_SNR1wedge30/wedge-mask.mrc'


import numpy as N
import tomominer.io.file as IF
v1 = IF.read_mrc_vol(v1f).astype(N.float)
m1 = IF.read_mrc_vol(m1f).astype(N.float)

if False:
    v2 = core.parse_mrc(v2f)
    m2 = core.parse_mrc(m2f)
else:
    v2 = v1
    m2 = m1

candidate_num = 200

angs = N.random.rand(candidate_num, 3) * N.pi * 2
locs_r = N.zeros([candidate_num, 3])


import tomominer.align.refine.gradient_refine as ARG
ARG.stocastic_parallel_search({'v1':v1, 'm1':m1, 'v2':v2, 'm2':m2, 'angs_org':angs, 'locs_r_org':locs_r, 'max_iter_num':100, 'stopping_tolerance':1e-3, 'min_dist_record_size':10, 'lambda':1, 'top_num':candidate_num/5})


'''


def fast_align_and_refine(op):
    if 'verbose' not in op:     op['verbose'] = False

    v1 = op['v1']
    m1 = op['m1'];      assert(m1.shape == v1.shape)
    v2 = op['v2'];    assert(v2.shape == v1.shape)
    m2 = op['m2'];    assert(m2.shape == v1.shape)

    if 'fast_max_l' in op:
        max_l = op['fast_max_l']
    else:
        max_l = 36    

    cs = core.combined_search(v1, m1, v2, m2, max_l) 
    cors = N.array([_[0] for _ in cs])        ;           assert  N.all(N.isfinite(cors))       ;       cors[cors > 1.0] = 1.0
    locs_r = N.array([_[1] for _ in cs])
    angs =   N.array([_[2] for _ in cs])

    cors_i = N.argsort(-cors)
    cors = cors[cors_i]
    locs_r = locs_r[cors_i]
    angs = angs[cors_i]

    if 'top_num' not in op:     op['top_num'] = len(cors)

    dists_sum = N.sqrt(2 - 2*cors)    # convert pearson correlation to euclidean distance

    candidate_num = len(dists_sum)
    
    sps_op = {'v1':v1, 'm1':m1, 'v2':v2, 'm2':m2, 'angs_org':angs, 'locs_r_org':locs_r, 'dists_sum':dists_sum, 'max_iter_num':(op['max_iter_num'] if 'max_iter_num' in op else 1000), 'stopping_tolerance':(op['stopping_tolerance'] if 'stopping_tolerance' in op else 1e-5), 'min_dist_record_size':(op['min_dist_record_size'] if 'min_dist_record_size' in op else 10), 'lambda':(op['lambda'] if 'lambda' in op else 1.0), 'top_num':(op['top_num'] if 'top_num' in op else N.inf), 'verbose':op['verbose']}

    if 'truth' in op:   sps_op['truth'] = copy.deepcopy(op['truth'])

    sps_re = stocastic_parallel_search(sps_op)

    #print 'WARNING: the following tests shows that there are some small difference between prediction and ground truth, dont know why'

    if False:
        # sort results according to dist_sum
        
        inds_i = sps_re['dist_sum'].argsort()

        re = {}
        re['dist_sum'] = sps_re['dist_sum'][inds_i]
        re['cor'] = sps_re['cor'][inds_i]
        re['ang'] = sps_re['ang'][inds_i]
        re['loc_r'] = sps_re['loc_r'][inds_i]
        
        return re
    else:
        return sps_re




'''
# test code


v1f = '/panfs/cmb-panasas2/mxu/proj/imaging/electron/frequent_structure/data/out/method/classification-system/pytom/thermosomes/particle/test_SNR1wedge30/test_SNR1wedge30_100.mrc'
m1f = '/panfs/cmb-panasas2/mxu/proj/imaging/electron/frequent_structure/data/out/method/classification-system/pytom/thermosomes/particle/test_SNR1wedge30/wedge-mask.mrc'
v2f = '/panfs/cmb-panasas2/mxu/proj/imaging/electron/frequent_structure/data/out/method/classification-system/pytom/thermosomes/particle/test_SNR1wedge30/test_SNR1wedge30_102.mrc'
m2f = '/panfs/cmb-panasas2/mxu/proj/imaging/electron/frequent_structure/data/out/method/classification-system/pytom/thermosomes/particle/test_SNR1wedge30/wedge-mask.mrc'

import numpy as N
import tomominer.io.file as IF
v1 = IF.read_mrc_vol(v1f).astype(N.float)
m1 = IF.read_mrc_vol(m1f).astype(N.float)

if False:
    v2 = IF.read_mrc_vol(v2f).astype(N.float)
    m2 = IF.read_mrc_vol(m2f).astype(N.float)
else:
    v2 = N.copy(v1)
    m2 = N.copy(m1)


import tomominer.align.refine.gradient_refine as ARG
faaf_re = ARG.fast_align_and_refine({'v1':v1, 'm1':m1, 'v2':v2, 'm2':m2, 'top_num':10, 'verbose':True})
print faaf_re



'''




'''
# test code

import tomominer.align.fast.util as AFU
tp = AFU.generate_toy_subtomogram_pair(dim_siz=32, snr=0.01)

import tomominer.align.refine.gradient_refine as ARG
faaf_re = ARG.fast_align_and_refine({'v1':tp['v1b'], 'm1':tp['m1'], 'v2':tp['v2b'], 'm2':tp['m2'], 'truth':{'ang':tp['ang'], 'loc_r':tp['loc']}, 'top_num':10, 'stopping_tolerance':1e-5, 'verbose':True})
print faaf_re

import tomominer.geometry.ang_loc as GA
print GA.angle_zyz_translation_difference(ang1=tp['ang'], loc1_r=tp['loc'], ang2=faaf_re['ang'], loc2_r=faaf_re['loc_r'])
print GA.angle_difference_zyz(ang1=tp['ang'], ang2=faaf_re['ang'])

import tomominer.io.file as IF
IF.put_mrc(tp['v1b'], '/tmp/mrc/v1b.mrc', overwrite=True)
IF.put_mrc(tp['v2b'], '/tmp/mrc/v2b.mrc', overwrite=True)


'''

