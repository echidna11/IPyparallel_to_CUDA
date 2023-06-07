
'''
use genetic algorithm to find a subset of subtomograms that maximizes ssnr based FSC
todo: when number of subtomograms is large, can first cluster (over partition) subtomograms, then evolve at cluster level
~/ln/tomominer/tomominer/pursuit/multi/pattern_generate/genetic_algorithm_ssnr_fsc.py

'''
'''
load subtomograms, and transform them
'''

import os, sys, copy, uuid, time
import numpy as N
import tomominer.statistics.ssnr as SS
import socket

def data_prepare(dj, op):
    print 'data_prepare()'

    if 'segmentation_tg_op' not in op:      op['segmentation_tg_op'] = None

    s_sum = None
    s_prod_sum = None
    s_mask_sum = None

    for i, d in enumerate(dj):
        print i, '      ', '\r',        ;       sys.stdout.flush()

        r = SS.var__local(self=None, data_json=[d], return_key=False, segmentation_tg_op=op['segmentation_tg_op'])

        if s_sum is None:
            siz = [len(dj)]
            siz.extend(r['sum'][0].shape)

            s_sum = N.zeros( siz, dtype=N.complex )
            s_prod_sum = N.zeros( siz, dtype=N.complex )
            s_mask_sum = N.zeros( siz )

        s_sum[i,:,:,:] = r['sum'][0]
        s_prod_sum[i,:,:,:] = r['prod_sum'][0]
        s_mask_sum[i,:,:,:] = r['mask_sum'][0]


    return {'sum':s_sum, 'prod_sum':s_prod_sum, 'mask_sum':s_mask_sum}




'''
Generic algorithm, iteratively evolve a pool of candidates
parameters:      stat: statistics prepared using data_prepare
'''
def ga(self=None, stat=None, initial_population=None, op=None):

    c_len = len(stat['sum'])            # the length of chromosome

    n = op['population_size']               # size of population

    if initial_population is not None:
        p = ga__init(p0=initial_population, size=n)
    else:
        assert op['init']['rand_threshold'] > 0
        assert op['init']['rand_threshold'] < 1

        p = (N.random.random( (op['population_size'], c_len) ) >= op['init']['rand_threshold']).astype(N.int)           # there are more ones when op['init']['rand_threshold'] is smaller

    if ('parallel' in op):
        if 'host' not in op['parallel']['redis']:       op['parallel']['redis']['host'] = socket.gethostname()
        stat_keys_key = ga__parallel__stat_distribute(self, stat, batch_id=str(uuid.uuid4()), op=op['parallel'])

    best_score = None

    best = None         # we keep a record of the best solutions occured so far

    for iter_i in range(op['max_iteration_num']):
        if 'parallel' in op:
            gep_op = copy.deepcopy(op['evaluate'])
            gep_op['parallel'] = copy.deepcopy(op['parallel'])

            e = ga_evaluate__parallel(self=self, p=p, stat_keys_key=stat_keys_key, op=gep_op)
        else:
            e = ga_evaluate(p=p, stat=stat, op=op['evaluate'])

        combined = ga__combine_population_ordered(c0=best, c1={'p':copy.deepcopy(p), 'e':copy.deepcopy(e)})         # put the best individuals together with current evaluated population

        p = ga_evolve(  p0=combined['p'], s=N.array( [_['score'] for _ in combined['e']] ), mr=op['mutation_rate'], n=n, sum_min=op['sum_min']  )

        best = {'p':combined['p'][:n,], 'e':combined['e'][:n]}

        # determine whether should stop
        max_i = N.argmax([_['score'] for _ in combined['e']])
        best_score_t = combined['e'][max_i]['score']

        if (best_score is None) or (best_score_t > best_score['score']):
            best_score = {'score':best_score_t, 'repeat':0}
        else:
            best_score['repeat'] += 1

        sel_nums = [N.sum(_) for _ in combined['p']]
        print '\r', '\t\t', 'iter_i %4d'%(iter_i,), 'score %3.4f'%(best_score_t, ), '        ', 'best_rep %3d'%(best_score['repeat'],), '        ', 'subtomogram num %6d'%(int(N.sum(combined['p'][max_i,:])),), '        ', 'min', N.min(sel_nums), 'max', N.max(sel_nums),            ;           sys.stdout.flush()

        if best_score['repeat'] >= op['best_score_max_repeat']:      break

        
    if ('parallel' in op):    ga__parallel__stat_cleanup(self=self, stat_keys_key=stat_keys_key)

    print

    return best



def ga__init(p0, size):

    p0 = copy.deepcopy(p0)

    p = N.zeros( (size, p0.shape[1]), dtype=N.int )
    for i in range(size):        p[i,:] = p0[i % len(p0), :]

    return N.array(p)
    


def ga_evaluate__single(l, stat, op):
    l = l.astype(bool)          # must use boolean for the indexing, not integer array!!!

    sum_v = N.sum(stat['sum'][l, :, :, :], axis=0)
    prod_sum = N.sum(stat['prod_sum'][l, :, :, :], axis=0)
    mask_sum = N.sum(stat['mask_sum'][l, :, :, :], axis=0)

    st = SS.ssnr__given_stat(sum_v=sum_v, prod_sum=prod_sum, mask_sum=mask_sum, op=op['ssnr'])
    st['ssnr_log_sum'] = N.sum(N.log(st['fsc']))
    st['fsc_sum'] = N.sum(st['fsc'])

    return st


def ga_evaluate__scoring(s, op):
    if op['method'] == 'fsc_sum':
        for ss in s:    ss['score'] = ss['fsc_sum']
    elif op['method'] == 'ssnr_log_sum':
        for ss in s:    ss['score'] = ss['ssnr_log_sum']
    else:
        raise Exception('method')


def ga_evaluate(p, stat, op):
    s = [None] * len(p)
    for i, l in enumerate(p):
        s[i] = ga_evaluate__single(l=l, stat=stat, op=op)

    ga_evaluate__scoring(s, op=op['scoring'])

    return s


'''
combine two populations and their scores, then order individuals according to scores
'''
def ga__combine_population_ordered(c0, c1):
    
    p = None
    e = None

    if c0 is not None:
        p = c0['p']
        e = c0['e']

        if c1 is not None:
            p = N.vstack(   (p, c1['p'])    )
            e.extend(c1['e'])

    elif c1 is not None:
        p = c1['p']
        e = c1['e']

    s = N.array([_['score'] for _ in e])
    i = N.argsort( -s )

    p = p[i,:]
    e = [e[_] for _ in i]

    return {'p':p, 'e':e}

 

def ga__parallel__stat_distribute(self, stat, batch_id, op):

    if ('redis' in op) and ('redis' not in dir(self)):
        from tomominer.io.redis_proxy import Redis
        self.redis = Redis(host=op['redis']['host'], port=op['redis']['port'], db=op['redis']['db'], expire_time=op['redis']['expire_time'])


    n_chunk = op['n_chunk']

    stat_keys = []        # a list of keys for stat data

    n = len(stat['sum'])
    ind = list(range(n))
    while ind:
        ind_t = ind[:n_chunk]
        stat_t = {'sum':{_:stat['sum'][_] for _ in ind_t}, 'prod_sum':{_:stat['prod_sum'][_] for _ in ind_t}, 'mask_sum':{_:stat['mask_sum'][_] for _ in ind_t}}

        stat_key_t = 'stat-key--' + str(uuid.uuid4())
        self.redis.set(stat_key_t, {'stat':stat_t, 'ind':ind_t})

        stat_keys.append(stat_key_t)

        ind = ind[n_chunk:]


    stat_keys_key = 'stat-keys--' + batch_id
    self.redis.set(stat_keys_key, stat_keys)

    return stat_keys_key


def ga__parallel__stat_cleanup(self, stat_keys_key):

    stat_keys = self.redis.get_del(stat_keys_key)

    for k in stat_keys:     self.redis.delete(k)


def ga_evaluate__parallel(self, p, stat_keys_key, op):
    
    n = p.shape[0]      # population size

    n_chunk = self.runner.estimate_chunk_size(n=n)
   
    ind = list(range(n))
    pl = p.tolist()

    tasks = []
    while ind:
        ind_t = ind[:n_chunk]
        pt = [pl[_] for _ in ind_t]

        population_key_t = 'population-key--' + str(uuid.uuid4())
        self.redis.set(population_key_t, {'ind':ind_t, 'p':pt})
    
        tasks.append(self.runner.task(module='tomominer.pursuit.multi.pattern_generate.genetic_algorithm_ssnr_fsc', method='ga_evaluate__parallel__local', kwargs={'population_key':population_key_t, 'stat_keys_key':stat_keys_key, 'op':op}))
        
        ind = ind[n_chunk:]


    re = [_.result for _ in self.runner.run__except(tasks)]

    s = [None] * p.shape[0]
    for rk in re:
        r = self.redis.get_del(rk)
        for ind_i, i in enumerate(r['ind']):
            s[i] = r['s'][ind_i]

    for ss in s:        assert ss is not None

    ga_evaluate__scoring(s, op=op['scoring'])

    return s



'''
given the following design, the code keep loading chunks of stat, which will lead to heavey io. Therefore it is good to use a very large population pool size (say 10000) to introduce more calculation
'''
def ga_evaluate__parallel__local(self, population_key, stat_keys_key, op):

    if 'verbose' not in op:     op['verbose'] = False

    func_start_time = time.time()

    if ('redis' in op['parallel']) and ('redis' not in dir(self)):        
        from tomominer.io.redis_proxy import Redis
        self.redis = Redis(host=op['parallel']['redis']['host'], port=op['parallel']['redis']['port'], db=op['parallel']['redis']['db'], expire_time=op['parallel']['redis']['expire_time'])

    
    io_cumulative_time = 0.0

    io_cur_time = time.time()
    pi_t = self.redis.get_del(population_key)
    io_cumulative_time += time.time() - io_cur_time

    ind = pi_t['ind']
    p = pi_t['p']
    del pi_t


    stat_keys = self.redis.get(stat_keys_key)

    sum_d = {}          # statistics of all chromosomes across all positions
    prod_sum_d = {}
    mask_sum_d = {}

    covered_positions = N.zeros(len(p[0]), dtype=N.int)


    for stat_key in stat_keys:
        #if self.work_queue.done_tasks_contains(self.task.task_id):            raise Exception('Duplicated task')           # dont use this inside this loop, becasue this check is very time consuming
        
        io_cur_time = time.time()
        stat_load = self.redis.get(stat_key)
        io_cumulative_time += time.time() - io_cur_time

        ind_pos = stat_load['ind']          # indics in different positions in a chromosome
        sum_t = stat_load['stat']['sum']
        prod_sum_t = stat_load['stat']['prod_sum']
        mask_sum_t = stat_load['stat']['mask_sum']

        covered_positions[ind_pos] = 1


        for i, l in enumerate(p):           # i is the relative index of individual, l is content in chromosome

            for ind_pos_t in ind_pos:
                if l[ind_pos_t] == 0:        continue            # if the corresponding subtomogram is excluded, ignore it

                if i not in sum_d:
                    sum_d[i] = N.zeros(sum_t[ind_pos_t].shape, dtype=N.complex)
                    prod_sum_d[i] = N.zeros(prod_sum_t[ind_pos_t].shape, dtype=N.complex)
                    mask_sum_d[i] = N.zeros(mask_sum_t[ind_pos_t].shape)

                sum_d[i] += sum_t[ind_pos_t]
                prod_sum_d[i] += prod_sum_t[ind_pos_t]
                mask_sum_d[i] += mask_sum_t[ind_pos_t]


    for _ in covered_positions:     assert _ == 1

    s = [None] * len(p)
    for i in sum_d:
        #if self.work_queue.done_tasks_contains(self.task.task_id):            raise Exception('Duplicated task')           # dont use this inside this loop, becasue this check is very time consuming

        st = SS.ssnr__given_stat(sum_v=sum_d[i], prod_sum=prod_sum_d[i], mask_sum=mask_sum_d[i], op=op['ssnr'])
        st['ssnr_log_sum'] = N.sum(N.log(st['fsc']))
        st['fsc_sum'] = N.sum(st['fsc'])
        s[i] = st

    for ss in s:        assert ss is not None

    io_cur_time = time.time()
    return_key = 'ga_evaluate__parallel__local' + str(uuid.uuid4())
    self.redis.set(return_key, {'ind':ind, 's':s})
    io_cumulative_time += time.time() - io_cur_time


    if op['verbose']:        print 'func time', time.time()-func_start_time, 'io time', io_cumulative_time

    #if self.work_queue.done_tasks_contains(self.task.task_id):            raise Exception('Duplicated task')

    return return_key




'''
paramters:  p0: given population,   s: score        mr: mutation rate       n: number of offsprings to generate     sum_min: the minimum number of candidates to be selected
'''

def ga_evolve(p0, s, mr, n=None, sum_min=None):
    assert sum_min is not None

    cdf = ga_evolve__make_ecdf(s)


    if n is None:       n = p0.shape[0]     # population pool size
    
    p = N.zeros( (n, p0.shape[1]) ) + N.nan

    c = 0
    while c < n:
        (i0, i1) = ga_evolve__select_pair(c=cdf)
        (l0, l1) = ga_evolve__crossover(p0[i0], p0[i1])

        pt = ga_evolve__mutate(l=l0, mr=mr)
        if pt.sum() < sum_min:  continue
        p[c,:] = pt
        c += 1

        if c >= n:  break

        pt = ga_evolve__mutate(l=l1, mr=mr)
        if pt.sum() < sum_min:  continue
        p[c,:] = pt
        c += 1

    return p.astype(N.int)



'''
get a emperical cumulative distribution function from scores
'''
def ga_evolve__make_ecdf(s):
    s = s - N.min(s)        # this is useful when the average s is much larger than the variance of s

    s_sum = N.sum(s)
    if s_sum > 0:
        s = s / float(s_sum)
    else:
        assert  N.isfinite(s_sum)
        s = N.ones(len(s)) * (1.0 / len(s))        # in such case, somehow all scores are equal, let's just uniform sampling

    c = 0
    cdf = N.zeros(len(s))
    for i in range(len(s)):
        c += s[i]
        cdf[i] = c

    return cdf



'''
sample from ecdf c
'''
def ga_evolve__sample_ecdf(c):
    return N.sum(c < N.random.uniform())




'''
select index of a pair of individuals according to some cumulative density function c
'''
def ga_evolve__select_pair(c, max_trial=1000):
    i0 = ga_evolve__sample_ecdf(c)

    c = 0
    while True:
        i1 = ga_evolve__sample_ecdf(c)
        if i0 != i1:    break

        c += 1
        assert c < max_trial

    return (i0, i1)
    



'''
crossover
'''

def ga_evolve__crossover(c0, c1):
    m = len(c0)

    i = N.random.randint(m)

    l0 = N.zeros(m, dtype=N.int) + N.nan
    l1 = N.zeros(m, dtype=N.int) + N.nan

    l0[:i] = c0[:i]
    l0[i:] = c1[i:]

    l1[:i] = c1[:i]
    l1[i:] = c0[i:]

    return (l0, l1)



def ga_evolve__mutate(l, mr):
    l = N.copy(l)

    for i in range(len(l)):
        if N.random.uniform() > mr:     continue

        if l[i] == 1:
            l[i] = 0
        elif l[i] == 0:
            l[i] = 1
        else:
            raise Exception()

    return l





'''

# test code using real data

data_file = '/home/rcf-47/mxu/ln/frequent_structure/data/out/analysis/shan/wang/150113/templates/matching/0000/pursuit/5-001-007/0000/recursive/clus/clus-003/data.json'

import json
with open(data_file) as f:   dj = json.load(f)


import tomominer.pursuit.multi.pattern_generate.genetic_algorithm_ssnr_fsc as PMPG
d = PMPG.data_prepare(dj[:400], op={'segmentation':None})


import cPickle as pickle
with open('/tmp/dp.pickle', 'wb') as f:     pickle.dump(d, f, protocol=-1)





import cPickle as pickle
with open('/tmp/dp.pickle', 'rb') as f:     d = pickle.load(f)

m = 200
d = {'sum':d['sum'][:m], 'prod_sum':d['prod_sum'][:m], 'mask_sum':d['mask_sum'][:m]}

import tomominer.pursuit.multi.pattern_generate.genetic_algorithm_ssnr_fsc as PMPG
PMPG.ga(stat=d, op={'population_size':50, 'max_iteration_num':1000, 'mutation_rate':0.01, 'best_score_max_repeat':50})




#--------------------------------------------------------------
# code for parallel computation


import cPickle as pickle
with open('/tmp/dp.pickle', 'rb') as f:     d = pickle.load(f)

#m = 1000
#d = {'sum':d['sum'][:m], 'prod_sum':d['prod_sum'][:m], 'mask_sum':d['mask_sum'][:m]}


from tomominer.parallel.queue_master import QueueMaster

import tomominer.common.obj as CO
from tomominer.io.cache import Cache
    
# just get a generic object that can contain a cache
self = CO.Object()

self.cache = Cache(tmp_dir='/home/rcf-47/mxu/tmp-panasas2')

self.runner = QueueMaster('localhost', 5011)

import tomominer.pursuit.multi.pattern_generate.genetic_algorithm_ssnr_fsc as PMPG
PMPG.ga(self=self, stat=d, op={'population_size':100, 'max_iteration_num':1000, 'mutation_rate':0.01, 'best_score_max_repeat':50, 'parallel':{'n_chunk':100}})


# todo: need to modify parameters according to /auto/cmb-08/fa/mxu/proj/imaging/electron/frequent_structure/data/out/analysis/jensen/chang_yiwei/140514/Bdellovibrio/classification/classification/pn/0005/pursuit/large-blob/0000/recursive/clus-000/avg/0002-parallel/average__op.json

'''


#------------------------------------------------------------------------------------------------------------

'''

# test code using simulation data. We want to know even if the subtomograms are from same structure of, say ribosome, and they are perfectly aligned, can GA tend to select all subtomograms to maximize SSNR?

#first generate simulation data, density map is randomly rotated, then convert to subtomogram, then rotate back. Save it to somewhere and generate data.json file. This is done in /home/rcf-47/mxu/ln/frequent_structure/data/out/method/ga-averaging/ribosome



wedge_angle = 30
snr = '0.01'

import os
data_json_in = '/home/rcf-47/mxu/ln/frequent_structure/data/out/method/ga-averaging/ribosome/subtomograms/%d/%s/data_config.json'%(wedge_angle,snr)

import json
with open(data_json_in) as f:   dj_in = json.load(f)

import tomominer.pursuit.multi.pattern_generate.genetic_algorithm_ssnr_fsc as PMPG
dj = PMPG.test_simulation_data_prepare(dj_in)



import tomominer.pursuit.multi.pattern_generate.genetic_algorithm_ssnr_fsc as PMPG
d = PMPG.data_prepare(dj, op={'segmentation':None})


import cPickle as pickle
with open('/tmp/dp-%d-%s.pickle'%(wedge_angle,snr), 'wb') as f:     pickle.dump(d, f, protocol=-1)





import cPickle as pickle
with open('/tmp/dp-30-0.01.pickle', 'rb') as f:     d = pickle.load(f)

import tomominer.pursuit.multi.pattern_generate.genetic_algorithm_ssnr_fsc as PMPG
r = PMPG.ga(        stat=d, op={        'population_size':50, 'max_iteration_num':1000, 'mutation_rate':0.01, 'sum_min':10, 'best_score_max_repeat':50, 'init':{'rand_threshold':0.5}, 'evaluate':{'ssnr':{'method':1, 'mask_sum_threshold bak':0}, 'scoring':{'method':'fsc_sum'}}      }         )




# observations: when wedge angle=30, and SNR=10, the method tend to pick up 60 out of 100 subtomograms. On the other hand, when SNR=0.1, the method tend to pick out 84 subtomograms. When SNR=0.01, the method tend to pick out 67 subtomograms. So the GA and SSNR-FSC score does not tend to select only half of particles due to any fault. YEAH !!!!


'''





'''
load data prepared by simulation.individual.subtomogram_generation
then rotate back subtomograms and corresponding masks so that they are aligned.
'''

def test_simulation_data_prepare(dj):

    import tomominer.geometry.ang_loc as AAL

    djn = []
    for d in dj:

        ang_inv, loc_inv = AAL.reverse_transform_ang_loc(N.array(d['model']['angle']), N.array(d['model']['loc']))

        djn.append(         {'subtomogram':d['subtomogram'], 'mask':d['mask'], 'angle':ang_inv.tolist(), 'loc':loc_inv.tolist()}          )


    return djn

