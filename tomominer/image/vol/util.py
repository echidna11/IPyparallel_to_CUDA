

'''
~/ln/tomominer/tomominer/image/vol/util.py
tomominer.image.vol.util
'''


import numpy as N




#-------------------------------------------------------------------------------------------------------------
# functions for geometrical manipulation

def roll(v, s0, s1, s2):
    if s0 != 0: 
        v0 = N.roll(v, s0, axis=0)
    else:
        v0 = v
    
    if s1 != 0:
        v1 = N.roll(v0, s1, axis=1)
    else:
        v1 = v0
    
    if s2 != 0: 
        v2 = N.roll(v1, s2, axis=2)    
    else:
        v2 = v1

    return v2



# rolling a volume for one step
def roll_one_step(v, s0, s1, s2):
    if s0 != 0: 
        #v0 = N.roll(v, s0, axis=0)
        v0 = N.zeros(v.shape, dtype=v.dtype)
        if s0 == 1:
            v0[1:v.shape[0], :, :] = v[0:(v.shape[0] - 1), :, :]
        elif s0 == -1:
            v0[0:(v.shape[0]-1), :, :] = v[1:v.shape[0], :, :]
        else:            assert False

    else:
        v0 = v
    
    if s1 != 0:
        #v1 = N.roll(v0, s1, axis=1)
        v1 = N.zeros(v.shape, dtype=v.dtype)
        if s1 == 1:
            v1[:, 1:v.shape[1], :] = v0[:, 0:(v.shape[1] - 1), :]
        elif s1 == -1:
            v1[:, 0:(v.shape[1]-1), :] = v0[:, 1:v.shape[1], :]
        else:       assert False

    else:
        v1 = v0
    
    if s2 != 0: 
        #v2 = N.roll(v1, s2, axis=2)    
        v2 = N.zeros(v.shape, dtype=v.dtype)
        if s2 == 1:
            v2[:, :, 1:v.shape[2]] = v1[:, :, 0:(v.shape[2] - 1)]
        elif s2 == -1:
            v2[:, :, 0:(v.shape[2]-1)] = v1[:, :, 1:v.shape[2]]
        else:       assert False

    else:
        v2 = v1

    return v2




# resize an volume(v) to given size, and keep image center same
import scipy.ndimage.interpolation as inter
def resize_center(v, s, cval=float('NaN')):

    vs = N.array(v.shape, dtype=N.float)
    
    v1 = inter.affine_transform(input=v, matrix=N.eye(v.ndim), offset=(vs-s)/2.0, output_shape=s, cval=cval )
    return v1


"""
# testing commands

import os,sys
sys.path.append(os.path.join( os.getenv('HOME'), 'ln/tomo/py' ))

import numpy as np
import tomominer.image.vol.util as gv
import tomominer.filter.gaussian as fg
v = fg.gauss_function([49,40,30], 5)
v1 = gv.resize_center(v, N.array(v.shape, dtype=N.float) / 2)
v2 = gv.resize_center(v, N.array(v.shape, dtype=N.float) * 2)


gv.dspcub(v)
gv.dspcub(v1)
gv.dspcub(v2)

import matplotlib.pyplot as plt
plt.show()


"""

# given a dictionary of volumes, find the largest, then generate a new set of volumes of same size as the largest multiplied by a factor
def resize_center_batch_dict(vs, cubic=True, enlarge_factor=None, cval=float('NaN')):
    siz = [N.array(vs[_].shape, dtype=int) for _ in vs]
    siz = N.array(siz)
    if cubic:
        siz = siz.max()
        siz = N.array([siz, siz, siz])
    else:
        siz = siz.max(axis=0)

    if enlarge_factor is not None:
        siz = N.ceil(siz * enlarge_factor).astype(int)

    vsn = {}
    for i in vs:
        vsn[i] = resize_center(vs[i], siz, cval)

    return vsn





def center_mass(v):
    size = v.shape
    x = N.mgrid[0:size[0], 0:size[1], 0:size[2]]

    v_sum = v.sum()

    cm = N.zeros(3)
    for dim_i in range(3):  cm[dim_i] = N.sum(x[dim_i] * v) / v_sum

    return cm



def fft_mid_co(siz):
    assert(all(N.mod(siz, 1) == 0))
    assert(all(N.array(siz) > 0))

    mid_co = N.zeros(len(siz), N.dtype("int64"))

    # according to following code that uses numpy.fft.fftshift()
    for i in range(len(mid_co)):
        m = siz[i]
        mid_co[i] = N.floor(m/2)

    return mid_co


'''
# verification

import os,sys
sys.path.append(os.path.join( os.getenv('HOME'), 'ln/tomominer/tomominer' ))


import numpy as N
import tomominer.image.vol.util as CV

for i in range(1000):
    siz = [N.random.randint(10, 70), N.random.randint(10, 70), N.random.randint(10, 70)]
    c = CV.fft_mid_co(siz)
    #print c
    v = N.zeros(siz)
    v[0,0,0] = 1
    v = N.fft.fftshift(v)
    #print N.where(v==1)
    if v[c[0], c[1], c[2]] == 1.0:      print i,




'''

#-----------------------------------------------------------------------------------

# roughly add a small vol to a big whole map, so that the center of vol is roughly centered at c
def add_to_whole_map(whole_map, vol, c=None):

    if c is None:   c = N.array(whole_map.shape) / 2

    c = N.round(c)

    siz = N.array(vol.shape)

    se = subvolume_center_start_end(c, map_siz=N.array(whole_map.shape), subvol_siz=siz)
    if se is None:       return None

    
    # we also handle NaN in the whole map, and replace them with valid values in local map (if any)
    local_map = whole_map[se[0,0]:se[0,1], se[1,0]:se[1,1], se[2,0]:se[2,1]]
    local_map[N.isnan(local_map)] = 0
    local_map += vol

    whole_map[se[0,0]:se[0,1], se[1,0]:se[1,1], se[2,0]:se[2,1]] = local_map

    return se


# roughly paste a small vol to a big whole map, so that the center of vol is roughly centered at c
def paste_to_whole_map(whole_map, vol, c=None):

    if c is None:   c = N.array(whole_map.shape) / 2

    c = N.round(c)

    siz = N.array(vol.shape)

    se = subvolume_center_start_end(c, map_siz=N.array(whole_map.shape), subvol_siz=siz)
    if se is None:       return None

    paste_to_whole_map__se(whole_map, vol, se)

    return se


# paste to a map given start and end coordinates
def paste_to_whole_map__se(whole_map, vol, se):
    whole_map[se[0,0]:se[0,1], se[1,0]:se[1,1], se[2,0]:se[2,1]] = vol


def cut_from_whole_map(whole_map, c, siz):

    se = subvolume_center_start_end(c, map_siz=N.array(whole_map.shape), subvol_siz=siz)
    return          cut_from_whole_map__se(whole_map, se)


# cut a map given start and end coordinates
def cut_from_whole_map__se(whole_map, se):
    if se is None:       return None
    return          whole_map[se[0,0]:se[0,1], se[1,0]:se[1,1], se[2,0]:se[2,1]]



# given a center c, get the relative start and end position of a subvolume with size subvol_siz
def subvolume_center_start_end(c, map_siz, subvol_siz):

    siz_h = N.ceil( subvol_siz / 2.0 )

    start = c - siz_h
    start = start.astype(int)
    end = start + subvol_siz
    end = end.astype(int)

    if any(start < 0):  return None
    if any(end >= map_siz):    return None

    se = N.zeros( (3,2), dtype=N.int )
    se[:,0] = start
    se[:,1] = end

    return se






#-----------------------------------------------------------------------------------
# grid of neighborhoods

# indicis of neighbor voxels
def neighbor_inds(i, dims):
    x = N.array(N.unravel_index(i, dims))

    for d0 in range(-1, 2):
        for d1 in range(-1, 2):
            for d2 in range(-1, 2):
                if all(N.array([d0, d1, d2]) == 0):    continue

                y = x + [d0, d1, d2]

                if any(y < 0):      continue
                if any(y >= dims):  continue
                    
                yield N.ravel_multi_index(y, dims)        


# return a grid of neighbors (in matrix form) for indexing
def neighbor_grid(dims):
    i = N.array(range(N.prod(dims)))
    x = N.array(N.unravel_index(i, dims))

    j = N.zeros([len(i), 26], dtype=int) - 1

    col_i = 0
    for d0 in range(-1, 2):
        for d1 in range(-1, 2):
            for d2 in range(-1, 2):
                if all(N.array([d0, d1, d2]) == 0):    continue

                yt = N.copy(x);     yt[0]+=d0;    yt[1]+=d1;    yt[2]+=d2

                # boundary check
                inds = (yt[0] >= 0) 
                for dim_i in range(1, 3):
                    inds = inds & (yt[dim_i] >= 0)
                for dim_i in range(3):
                    inds = inds & (yt[dim_i] < dims[dim_i])    
                
                j[inds,col_i] = N.ravel_multi_index(yt[:,inds], dims)

                col_i += 1

    assert (col_i == 26), col_i

    j_filter = [None] * len(j)          # only a simple way to create an list with a fixed length
    for j_i in range(len(j)):       j_filter[j_i] = j[j_i, (j[j_i,:]>=0)]

    return j_filter
                

#-------------------------------------------------------------------
# grid functions

def grid_displacement_to_center(size, mid_co=None):

    size = N.array(size, dtype=N.float)
    assert size.ndim == 1

    if mid_co is None:            mid_co = (N.array(size) - 1) / 2          # IMPORTANT: following python convension, in index starts from 0 to size-1!!! So (siz-1)/2 is real symmetry center of the volume
 
    if size.size == 3:
        # construct a gauss function whose center is at center of volume
        grid = N.mgrid[0:size[0], 0:size[1], 0:size[2]]

        for dim in range(3):
            grid[dim, :, :, :] -= mid_co[dim]

    elif size.size == 2:
        # construct a gauss function whose center is at center of volume
        grid = N.mgrid[0:size[0], 0:size[1]]

        for dim in range(2):
            grid[dim, :, :] -= mid_co[dim]

    else:
        assert False

    return grid

def grid_distance_sq_to_center(grid):
    dist_sq = N.zeros(grid.shape[1:])
    if grid.ndim == 4:
        for dim in range(3):
            dist_sq += N.squeeze(grid[dim, :, :, :]) ** 2
    elif grid.ndim == 3:
        for dim in range(2):
            dist_sq += N.squeeze(grid[dim, :, :]) ** 2
    else:
        assert False

    return dist_sq

def grid_distance_to_center(grid):
    dist_sq = grid_distance_sq_to_center(grid)
    return N.sqrt(dist_sq)

#--------------------------------------------------------------------
# functions for visualization


# convert a 3D cube to a 2D image of slices
def cub_img(v, view_dir=2):
    if view_dir == 0:
        vt = N.transpose(v, [1,2,0])
    elif view_dir == 1:
        vt = N.transpose(v, [2,0,1])
    elif view_dir == 2:
        vt = v
    
    row_num = vt.shape[0] + 1
    col_num = vt.shape[1] + 1
    slide_num = vt.shape[2]
    disp_len = int( N.ceil(N.sqrt(slide_num)) )
    
    slide_count = 0
    im = N.zeros( (row_num*disp_len, col_num*disp_len) ) + float('nan')
    for i in range(disp_len):
        for j in range(disp_len):
            im[(i*row_num) : ((i+1)*row_num-1),  (j*col_num) : ((j+1)*col_num-1)] = vt[:,:, slide_count]
            slide_count += 1
            
            if (slide_count >= slide_num):
                break
            
        
        if (slide_count >= slide_num):
            break
   
    
    im_v = im[N.isfinite(im)]

    if im_v.max() > im_v.min(): 
        im = (im - im_v.min()) / (im_v.max() - im_v.min())

    return {'im':im, 'vt':vt}

# display an image
def dsp_img(v, new_figure=True):

    import matplotlib.pyplot as plt

    if new_figure:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        ax = plt


    import matplotlib.cm as cm
    
    ax_u = ax.imshow(  v, cmap = cm.Greys_r )
    ax.axis('off') # clear x- and y-axes

    plt.pause(0.001)        # calling pause will display the figure without blocking the program, see segmentation.active_contour.morphsnakes.evolve_visual



def dsp_cub(v, view_dir=2, new_figure=True):

    dsp_img(cub_img(v=v, view_dir=view_dir)['im'])



'''

# test code

import matplotlib
matplotlib.use('Qt4Agg')

import tomominer.model.util as MU
t = MU.generate_toy_model()

import tomominer.image.vol.util as IVU
IVU.dsp_cub(t)

'''


'''
displace three orthogonal slices of a volume
'''
def dsp_orthogonal_slices__matavi(v, c=None, vmin=None, vmax=None):

    if c == None:        c = N.array(v.shape) / 2.0
    
    c = tuple(c)

    from mayavi import mlab
    src = mlab.pipeline.scalar_field(v)

    if vmin is None:     vmin = v.min()
    if vmax is None:    vmax = v.max()

    cut_plane0 = mlab.pipeline.scalar_cut_plane(src, plane_orientation='x_axes', colormap='gray', vmin=vmin, vmax=vmax)
    cut_plane0.implicit_plane.origin = c
    cut_plane0.implicit_plane.widget.enabled = False

    cut_plane1 = mlab.pipeline.scalar_cut_plane(src, plane_orientation='y_axes', colormap='gray', vmin=vmin, vmax=vmax)
    cut_plane1.implicit_plane.origin = c
    cut_plane1.implicit_plane.widget.enabled = False

    cut_plane2 = mlab.pipeline.scalar_cut_plane(src, plane_orientation='z_axes', colormap='gray', vmin=vmin, vmax=vmax)
    cut_plane2.implicit_plane.origin = c
    cut_plane2.implicit_plane.widget.enabled = False


