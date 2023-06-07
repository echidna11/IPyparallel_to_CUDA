'''
Code for running pursuit code for individually simulted subtomograms
'''
import ipyparallelUtil as CU
import os, sys, json, copy, warnings
from tomominer.io.cache import Cache
import tomominer.common.obj as CO
import ipyparallel

def run_single(op):
    warnings.filterwarnings('error')
    self = CO.Object()
    self.cache = Cache(tmp_dir = op['options']['tmp_dir'])
    with open(op['data_file']) as f:
        dj = json.load(f)
    for d in dj:
        if not os.path.isabs(d['subtomogram']):
            d["subtomogram"] = os.path.abspath(os.path.join(os.path.dirname(op['data_file']), d['subtomogram']))
        if not os.path.isabs(d['mask']):
            d['mask'] = os.path.abspath(os.path.join(os.path.dirname(op['data_file']), d['mask']))
    #from tomominer.pursuit.multi.main_ipyparallel import pursuit
    import main_ipyparallel
    fs = main_ipyparallel.pursuit(self=self, op=op, data_json=dj)['file_stat_file']
    return {'file_stat': fs}

def main():
    rc = ipyparallel.Client()
    print(len(rc.ids))
    if len(rc.ids)==0:
        return
    with open('pursuit-combine-op.json') as f:
        op = json.load(f)
    out_dir = os.path.abspath(op["out_dir"])
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    results = []
    results.append(run_single(op))

if __name__ == '__main__':
    main()
