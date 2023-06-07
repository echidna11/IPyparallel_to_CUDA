# given an alignment against template file, plot all alignments aligned template-subtomogram pairs
# can order these pairs according to alignment scores

import sys
import json
import numpy as NP

import tomo
import general_util.vol as GV
import tomominer.io.file as FIO
import tomominer.filter.gaussian as FG

import matplotlib.pyplot as PLT
import matplotlib.cm as CM
import pylab



if __name__ == '__main__' :
    
    gauss_sigma = float(sys.argv[1])
    data_file = sys.argv[2]

    with open(data_file) as f:      data_json = json.load(f)

    # sort records according to alignment scores
    scores = NP.zeros(len(data_json))
    for i,rec in enumerate(data_json):      scores[i] = rec['score']

    si = NP.argsort(-scores)
    data_json_s = [data_json[i] for i in si]


    # export plots
    for sub_i, rec in enumerate(data_json_s):
        sys.stdout.write('%05d   %f   \r'%(sub_i, rec['score']));      sys.stdout.flush()

        tem = FIO.get_mrc(str(rec['template']))
        sub = FIO.get_mrc(str(rec['subtomogram']))

        sub_r = tomo.rotate_vol_pad_mean_py(sub, NP.array( rec['angle'], dtype=NP.float ), NP.array( rec['loc'], dtype=NP.float ))

        sub_rg = FG.smooth(v=sub_r, sigma=gauss_sigma)


        # plotting and storing figures
        f, sp = PLT.subplots(1,2)

        sp[0].imshow( GV.cub_img(v=tem)['im'] , cmap = CM.Greys_r )
        sp[1].imshow( GV.cub_img(v=sub_rg)['im'] , cmap = CM.Greys_r )

        for i in range(sp.size):
            sp[i].axis('off') # clear x- and y-axes

        PLT.draw()

        fig_file = '%05d.png'%(sub_i)
        pylab.savefig( fig_file )
        
        PLT.close("all")

      

