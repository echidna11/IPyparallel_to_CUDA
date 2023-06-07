#!/usr/bin/env python


'''
main routine
after running pursuit.multi.main.pursuit() to discover multiple patterns
then partition subtomograms according to best alignment
then, for each partition, call pursuit.multi.main.pursuit() to discover single pattern
'''

import os, sys, json, copy
import warnings
from multiprocessing.pool import Pool
from tomominer.parallel.queue_master import QueueMaster

import tomominer.common.obj as CO
from tomominer.io.cache import Cache





def main():

    warnings.filterwarnings('error')

    with open('./pursuit-op.json') as f:       op = json.load(f)

    op['template']['match']['priority'] = 10000         # just make the template matching part faster, because it does not cost much computation

     # first load file stat obtained from multi pursuit
    with open(os.path.join(op['out_dir'], 'file_stat.json')) as f:       file_stat = json.load(f)
    file_stat['passes'] = {int(_):file_stat['passes'][_] for _ in file_stat['passes']}

   
    # just get a generic object that can contain a cache
    self = CO.Object()

    self.pool = None
  
    self.cache = Cache(tmp_dir=op['options']['tmp_dir'])

    self.runner = QueueMaster(op['options']['network']['qhost'], op['options']['network']['qport'])

    if op['single_pursuit']['filter']['second_largest_cut']:
        import tomominer.pursuit.multi.recursive.filtering.second_largest_cut as PMRFS
        t = PMRFS.load_data(fs=file_stat)
        data_json = PMRFS.do_filter(al=t['al'], dj=t['dj'])
        del t
    else:
        with open(file_stat['passes'][file_stat['pass_i_current']]['data_json_file']) as f:      data_json = json.load(f)

    file_stat = pursuit(self=self, op=op, data_json=data_json, file_stat=file_stat)




def pursuit(self, data_json, file_stat, op):

    op['out_dir'] = os.path.abspath(op['out_dir'])

    # partition data
    single_out_dir = os.path.join(file_stat['passes'][file_stat['pass_i_current']]['pass_dir'], 'recursive', 'pursuit', 'single')

    if not os.path.isdir(single_out_dir):        os.makedirs(single_out_dir)
    file_stat['single_out_dir'] = single_out_dir

    psm_op_s = single_pursuit_data_op_prepare(fs=file_stat, op=op, out_dir=os.path.join(single_out_dir, 'sets'))

    file_stat['passes'][file_stat['pass_i_current']]['single'] = {}

    # perform single pursuit for each data
    import tomominer.pursuit.single.main as PSM
    for i in psm_op_s:
        with open(psm_op_s[i]['data_json_file']) as f:     dj_t = json.load(f)
        psm_file_stat = PSM.pursuit(self=self, op=psm_op_s[i], data_json=dj_t)
        file_stat['passes'][file_stat['pass_i_current']]['single'][i] = psm_file_stat['file_stat_file']

    with open(file_stat['file_stat_file'], 'w') as f:       json.dump(file_stat, f, indent=2)
    return file_stat




'''
given result of multi pursuit, we prepare data for single pursuit
IMPORTANT: do not overwrite existing files!!
'''
def single_pursuit_data_op_prepare(fs, op, out_dir):

    ops = {}

    fsp = fs['passes'][fs['pass_i_current']]
    with open(fsp['data_json_file']) as f:       dj = json.load(f)


    # first collect ids of matched templates
    templates = {}
    for _ in dj:
        if 'template' not in _:     continue
        if _['template']['id'] in templates:        continue
        templates[_['template']['id']] = _['template']

    for id_t in templates:
        out_dir_t = os.path.join(out_dir, '%d'%(id_t,))
        if not os.path.isdir(out_dir_t):        os.makedirs(out_dir_t)

        data_json_t_file = os.path.join(out_dir_t, 'data.json')
        if not os.path.isfile(data_json_t_file):
            with open(data_json_t_file, 'w') as f:      json.dump([_ for _ in dj if ('template' in _) and (_['template']['id'] == id_t)], f, indent=2)

        op_t_file = os.path.join(out_dir_t, 'pursuit-op.json')
        if os.path.isfile(op_t_file):
            with open(op_t_file):        op_t = json.load(f)
            assert      op_t['data_json_file'] == data_json_t_file
            assert      op_t['out_dir'] == os.path.join(out_dir_t, 'out')
        else:
            op_t = copy.deepcopy(op)
            op_t['data_json_file'] = data_json_t_file
            op_t['out_dir'] = os.path.join(out_dir_t, 'out')

        op_t['initial_template'] = templates[id_t]

        ops[id_t] = op_t

    return ops





if __name__ == '__main__':
    main()

