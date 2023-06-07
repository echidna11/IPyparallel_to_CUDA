
'''

just some simple tests of alignments


~/ln/tomominer/tomominer/pursuit/multi/tests/align/test_0.py
'''




v = '/auto/cmb-08/fa/mxu/proj/imaging/electron/frequent_structure/data/out/method/imputation-classification/tomogram-subtomogram/0003/subtomograms/5000/600-600-200/1/30/10/peak-dog/30/subtomograms/1/6/001651.mrc'
vm = '/auto/cmb-08/fa/mxu/proj/imaging/electron/frequent_structure/data/out/method/imputation-classification/tomogram-subtomogram/0003/subtomograms/5000/600-600-200/1/30/10/peak-dog/30/subtomograms/wedge-mask.mrc'

t = '/auto/cmb-08/fa/mxu/proj/imaging/electron/frequent_structure/data/out/method/imputation-classification/tomogram-subtomogram/0003/classification/0020/out/pass/004/common_frame/clus_vol_avg_007.mrc'
tm = '/auto/cmb-08/fa/mxu/proj/imaging/electron/frequent_structure/data/out/method/imputation-classification/tomogram-subtomogram/0003/classification/0020/out/pass/004/common_frame/clus_mask_avg_007.mrc'


import tomominer.io.file as IF
import tomominer.pursuit.multi.util as PMU

v = IF.get_mrc(v)
vm = IF.get_mrc(vm)

a = PMU.align_to_templates__pair_align(c=0, t_key={'subtomogram':t, 'mask':tm}, v=v, vm=vm, align_op={'L':36, 'with_missing_wedge':True})

print a


import numpy as N
import tomominer.core as core
core.write_mrc(v, '/tmp/v.mrc')


'''
2015-05-13 11:29:34,536 hpc2212          12105816         queue_worker WARNING  alignment failed: rec {u'loc': [-0.0, -1.0, -0.0], u'angle': [-1.3089969389957472, -2.443460952792061, -1.3089969389957472], u'mask': u'/auto/cmb-08/fa/mxu/proj/imaging/electron/frequent_structure/data/out/method/imputation-classification/tomogram-subtomogram/0003/subtomograms/5000/600-600-200/1/30/10/peak-dog/30/subtomograms/wedge-mask.mrc', u'score': 0.6337670858634568, u'template': {u'segmentation': u'/auto/cmb-08/fa/mxu/proj/imaging/electron/frequent_structure/data/out/method/imputation-classification/tomogram-subtomogram/0003/classification/0020/out/pass/003/common_frame/clus_vol_seg_phi_004.mrc', u'mask': u'/auto/cmb-08/fa/mxu/proj/imaging/electron/frequent_structure/data/out/method/imputation-classification/tomogram-subtomogram/0003/classification/0020/out/pass/003/common_frame/clus_mask_avg_004.mrc', u'id': 4, u'cluster': 23, u'subtomogram': u'/auto/cmb-08/fa/mxu/proj/imaging/electron/frequent_structure/data/out/method/imputation-classification/tomogram-subtomogram/0003/classification/0020/out/pass/003/common_frame/clus_vol_avg_004.mrc', u'pass_i': 1}, u'subtomogram': u'/auto/cmb-08/fa/mxu/proj/imaging/electron/frequent_structure/data/out/method/imputation-classification/tomogram-subtomogram/0003/subtomograms/5000/600-600-200/1/30/10/peak-dog/30/subtomograms/1/6/001651.mrc'}, template {'segmentation': u'/auto/cmb-08/fa/mxu/proj/imaging/electron/frequent_structure/data/out/method/imputation-classification/tomogram-subtomogram/0003/classification/0020/out/pass/004/common_frame/clus_vol_seg_phi_007.mrc', 'mask': u'/auto/cmb-08/fa/mxu/proj/imaging/electron/frequent_structure/data/out/method/imputation-classification/tomogram-subtomogram/0003/classification/0020/out/pass/004/common_frame/clus_mask_avg_007.mrc', 'id': 7, 'cluster': 22, 'subtomogram': u'/auto/cmb-08/fa/mxu/proj/imaging/electron/frequent_structure/data/out/method/imputation-classification/tomogram-subtomogram/0003/classification/0020/out/pass/004/common_frame/clus_vol_avg_007.mrc', 'pass_i': 1}, error None 

'''
