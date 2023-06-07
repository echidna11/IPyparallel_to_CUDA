

# functions to test different segmentation methods

import time
import sys
import numpy as np

def save_figure(fname, a_re):

    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    f, sp = plt.subplots(2,3)

    import general_util.vol as gv

    sp[0,0].imshow( gv.cub_img(v=a_re['vg'])['im'] , cmap = cm.Greys_r )
    sp[0,1].imshow( gv.cub_img(v=a_re['t_seg_r'])['im'] , cmap = cm.Greys_r )
    sp[0,2].imshow( gv.cub_img(v=a_re['sc_re']['stru_lbl'])['im'] , cmap = cm.Greys_r )
    sp[1,0].imshow( gv.cub_img(v=a_re['sc_re']['conn_rgn_lbl'])['im'] , cmap = cm.Greys_r )
    sp[1,1].imshow( gv.cub_img(v=a_re['sc_re']['mask'])['im'] , cmap = cm.Greys_r )
    sp[1,2].imshow( gv.cub_img(v=a_re['v_seg'])['im'] , cmap = cm.Greys_r )

    for i in range(sp.shape[0]):
        for j in range(sp.shape[1]):
            sp[i,j].axis('off') # clear x- and y-axes

    plt.draw()

    import pylab
    pylab.savefig(fname)
    
    plt.close("all")


if __name__ == '__main__':



    import os


    code_dir = os.path.abspath( os.path.join( os.path.dirname( os.path.abspath(__file__) ), '..' ) )

    import sys
    sys.path.append(code_dir)

    #import tomominer.plot     # for initializing remote display

    
    workspace_file = '/home/rcf-47/mxu/tmp1/em/beck/u2os/classification/subtomograms/1/_pursuit/pass_004/dump.pickle'

    import pickle
    with open(workspace_file) as f:      ws = pickle.load(f)

    #print ws.keys()


    vmal = ws['data']
    
    tk = ws['gavg_re']['vol_avg_out_key']
    tmk = ws['gavg_re']['mask_avg_out_key']
    
    import tomominer.io.vol as iv

    t = iv.get_mrc(tk)      # template
    tm = iv.get_mrc(tmk)    # mask

    gauss_sigma = 2

    import tomominer.filter.gaussian as fg
    tg = np.array( fg.smooth(v=t, sigma=gauss_sigma), dtype=np.float, order='F' )

    import general_util.vol as gv
    #gv.dspcub(tg)


    import tomominer.segmentation.thresholding as segt
    cutoff_num = 100
    cutoffs = np.linspace(0, tg.max(), cutoff_num)
    cutoffs = cutoffs[1:];      cutoffs = cutoffs[:len(cutoffs)-1]
    mmb = segt.ms_mo_best( segt.ms_cutoffs(v=tg, cutoffs=cutoffs) )
    t_seg =     np.array( (tg > mmb['best']['cutoff']), dtype=np.float, order='F' ) 
    #gv.dspcub(t_seg)

    #import matplotlib.pyplot as plt;    plt.show()

    import tomominer.align.util as au
    out_dir = os.path.join(os.getenv('HOME'), 'tmp1/em/target-complex-extraction/template-guided')
    for i, (vk, vmk, ang, loc) in enumerate(vmal):
        v = iv.get_mrc(vk)      # subtomogram
        vm = iv.get_mrc(vmk)    # mask
        vg = np.array( fg.smooth(v=v, sigma=gauss_sigma), dtype=np.float, order='F' )

        L = 36

        if False:
            import pickle;  
            with open('/tmp/v', 'wb') as f:     pickle.dump({'t':t, 'tg':tg, 't_seg':t_seg, 'tm':tm, 'v':v, 'vg':vg, 'vm':vm, 'L':L}, f)

        start_time = time.time()

        a_re = au.segment_and_align_one_vol_against_one_template(t=t, tg=tg, t_seg=t_seg, tm=tm, v=v, vg=vg, vm=vm, L=L)

        sys.stdout.write('%d,%d, %.2f   %.2f\t'%(i, a_re['sc_re']['conn_rgn_lbl'].max(), a_re['sc_re']['cutoff'] / a_re['sc_re']['max_intensity'], time.time() - start_time));     sys.stdout.flush()       # print out the number of connected regions detected, and time cost

        a_re['vg'] = vg
        
        save_figure(os.path.join(out_dir, '%05d.png'%(i)), a_re)


"""

import pickle
with open('/tmp/v') as f:   o=pickle.load(f)

import tomominer.align.util as au
a_re = au.segment_and_align_one_vol_against_one_template(t=o['t'], tg=o['tg'], t_seg=o['t_seg'], tm=o['tm'], v=o['v'], vg=o['vg'], vm=o['vm'], L=o['L'])

import os
import tomominer.segmentation.test as st
out_dir = os.path.join(os.getenv('HOME'), 'tmp2')
i=0
st.save_figure(os.path.join(out_dir, '%05d.png'%(i)), a_re)

"""


