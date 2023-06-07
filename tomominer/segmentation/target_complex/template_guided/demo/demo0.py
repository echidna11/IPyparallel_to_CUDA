#!/usr/bin/env python



'''

demostration of template guided target complex segmentation

~/ln/tomominer/tomominer/segmentation/target_complex/template_guided/demo/demo0.py

'''

import os, json, copy
import numpy as N


import tomominer.model.util as MU
import tomominer.image.optics.ctf as IOC
import tomominer.image.vol.util as IVU
import tomominer.filter.gaussian as FG


'''
generate a template of given size s
'''
def toy_template(s):

    v = N.zeros((s,s,s))

    s2 = s/2.0
    s4 = s/4.0
    s6 = s/6.0
    s8 = s/8.0

    v += MU.sphere_mask(shape=v.shape, center=(s2-s6, s2, s2), radius=s8)
    v += MU.sphere_mask(shape=v.shape, center=(s2+s6, s2, s2), radius=s8)

    
    return v




'''
generate a toy subtomogram of a give size s
'''
def toy_subtomogram(s):
    
    v = N.zeros((s,s,s))

    s2 = s/2.0
    s4 = s/4.0
    s6 = s/6.0
    s8 = s/8.0
    s12 = s/12.0

    v += MU.rectangle_mask(shape=v.shape, p0=(s2-s6-s8,s2-s8,s2-s8), p1=(s2-s6+s8, s2+s8, s2+s8))
    v += MU.rectangle_mask(shape=v.shape, p0=(s2+s6-s8,s2-s8,s2-s8), p1=(s2+s6+s8, s2+s8, s2+s8))

    v += MU.rectangle_mask(shape=v.shape, p0=(s2-s6,s2+s4+s12-s12,s2-s12), p1=(s2+s6, s2+s4+s12+s12, s2+s12)) * 0.8

    return v


def convolve_ctf(v, ctf_op):
    ctf = IOC.create(Dz=ctf_op['Dz'], size=v.shape, pix_size=ctf_op['pix_size'], voltage=ctf_op['voltage'], Cs=ctf_op['Cs'], sigma=(None if 'sigma' not in ctf_op else ctf_op['sigma']))['ctf']
    
    from numpy.fft import fftn, ifftn, fftshift, ifftshift

    return N.real(ifftn(ifftshift( fftshift(fftn(v)) * ctf )))


def main():
    with open('segmentation_target_complex_template_guided_demo_demo0__op.json') as f:     op = json.load(f)

    import matplotlib
    matplotlib.use('Qt4Agg')


    vt = toy_template(op['model']['size'])
    if not op['density_positive']:      vt = -vt
    vt = FG.smooth(v=vt, sigma=op['model']['smooth_sigma'])
    vt = convolve_ctf(vt, op['ctf'])
    #IVU.dsp_cub(vt)

    vs = toy_subtomogram(op['model']['size'])
    if not op['density_positive']:      vs = -vs
    vs = FG.smooth(v=vs, sigma=op['model']['smooth_sigma'])
    if N.isfinite(op['model']['SNR']):     vs += N.random.normal(scale=N.sqrt(vs.std() / op['model']['SNR']), size=vs.shape)
    vs = convolve_ctf(vs, op['ctf'])
    #IVU.dsp_cub(vs)

    import tomominer.pursuit.multi.util as PMU
    t_seg_op = copy.deepcopy(op['segmentation'])
    t_seg_op['density_positive'] = op['density_positive']
    t_seg_op['normalize_and_take_abs'] = op['template']['segmentation']['normalize_and_take_abs']
    t_phi = PMU.template_segmentation__single_vol(v=vt, op=t_seg_op)
    t_phi_mask = (t_phi>0.5)    ;       assert  t_phi_mask.sum() > 0        ;       assert N.logical_not(t_phi_mask).sum() > 0


    s_seg_op = copy.deepcopy(op['segmentation'])
    s_seg_op['density_positive'] = op['density_positive']
    s_seg_op['gaussian_smooth_sigma'] = op['template']['guided_segmentation']['gaussian_smooth_sigma']
    s_seg_op['phi_propotion_cutoff'] = op['template']['guided_segmentation']['phi_propotion_cutoff']
    vs_s = PMU.template_guided_segmentation(v=vs, m=t_phi_mask, op=s_seg_op)
    vs_s[N.logical_not(N.isfinite(vs_s))] = vs_s[N.isfinite(vs_s)].mean()
    

    #-------------------------------------------------
    # draw figures

    if not os.path.isdir(op['out_dir']):        os.makedirs(op['out_dir'])

    import matplotlib.pyplot as PLT
    import matplotlib.cm as MM

    fig = PLT.figure()
    PLT.imshow(vt[:, :, vt.shape[2]/2], cmap = MM.Greys_r )
    PLT.axis('off')
    PLT.draw()
    #PLT.show()
    PLT.savefig(os.path.join(op['out_dir'], 'vt.png'), bbox_inches='tight')
    PLT.close(fig)



    # plot the subtomogram, overlay with segment region
    fig = PLT.figure()
    img_t = vs[:, :, vs.shape[2]/2]       ;   img_t -= img_t.min()        ;       img_t /= img_t.max()
    img_tm = N.array(img_t)      ;           img_tm[t_phi_mask[:,:,t_phi_mask.shape[2]/2]] = 0.0
    img = N.zeros((vs.shape[0], vs.shape[1], 3))
    img[:,:,0] = img_t     ;           img[:,:,1] = img_tm      ;           img[:,:,2] = img_tm
    PLT.imshow(img)
    PLT.axis('off')
    PLT.draw()
    #PLT.show()
    PLT.savefig(os.path.join(op['out_dir'], 'vs.png'), bbox_inches='tight')
    PLT.close(fig)


    fig = PLT.figure()
    PLT.imshow(vs_s[:, :, vs_s.shape[2]/2], cmap = MM.Greys_r )
    PLT.axis('off')
    PLT.draw()
    #PLT.show()
    PLT.savefig(os.path.join(op['out_dir'], 'vs-s.png'), bbox_inches='tight')
    PLT.close(fig)



if __name__ == '__main__':
    main()



