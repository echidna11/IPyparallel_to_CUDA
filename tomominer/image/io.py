
# functions for loading and saving images

import numpy as N

def save_png(m, name, verbose=False):

    if verbose:
        print 'save_png()'
        print 'unique values', sorted(set(m.flatten()))

    m = N.array(m, dtype=N.float)
    mv = m[N.isfinite(m)]
    m = (m - mv.min()) / (mv.max() - mv.min())

    m = N.ceil(m * 65534)
    m = N.array(m, dtype=N.uint16)


    import png          # in pypng package
    png.from_array(m, 'L').save(name)




def save_image_matplotlib(m, out_file, vmin=None, vmax=None):
    import matplotlib.pyplot as PLT
    import matplotlib.cm as CM

    if vmin is None:        vmin = m.min()
    if vmax is None:        vmax = m.max()

    ax_u = PLT.imshow(  m, cmap = CM.Greys_r, vmin=vmin, vmax=vmax)
    PLT.axis('off')
    PLT.draw()

    PLT.savefig(out_file, bbox_inches='tight')
    PLT.close("all")


