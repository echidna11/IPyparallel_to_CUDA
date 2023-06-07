import numpy as N
import scipy.spatial.distance as SSD
import tomominer.io.file as IF
import math, os

def gaussian_value(x, sigma):
    tmp = (x*x)/(2*sigma*sigma)
    tmp = math.exp(-tmp)
    return tmp

def create_subtomogram(size, cords, rad, file_name, rescale=20):
    cords = cords/rescale
    rad = rad/rescale
    rad_siz = [N.int(10*N.ceil(_)) for _ in rad]
    size = N.int(N.ceil(size/rescale))
    m = N.zeros((2*size, 2*size, 2*size))
    for bead_i in range(len(cords)):
        c = N.array(N.ceil(cords[bead_i]), dtype=N.int).tolist()
        tmp_sphere = N.zeros((rad_siz[bead_i], rad_siz[bead_i], rad_siz[bead_i]))
        center_cord = tmp_sphere.shape[0]/2
        for i in range(tmp_sphere.shape[0]):
            for j in range(tmp_sphere.shape[1]):
                for k in range(tmp_sphere.shape[2]):
                    tmp_dis = SSD.euclidean([center_cord, center_cord, center_cord], [i+0.5, j+0.5, k+0.5])
                    if tmp_dis <= 1.5*rad[bead_i]:
                        tmp_sphere[i, j, k] = 1.0
                    elif tmp_dis > 3*rad[bead_i]:
                        tmp_sphere[i, j, k] = 0.0
                    else:
                        sigma_i = rad[bead_i]/3.0
                        t1 = gaussian_value((tmp_dis - 1.5*rad[bead_i]), sigma_i)
                        tmp_sphere[i, j, k] = t1
        
        box_size = tmp_sphere.shape[0]
        tmp = box_size/2
        t1 = range(tmp)
        t1.reverse()
        t1 = [-_ for _ in t1]
        t1.extend(range(tmp+1)[1:])
        assert len(t1) == N.int(box_size)
        cx = [c[0]+_+size for _ in t1]
        cy = [c[1]+_+size for _ in t1]
        cz = [c[2]+_+size for _ in t1]
        cx = [_ for _ in cx if _>=0 and _<m.shape[0]]
        cy = [_ for _ in cy if _>=0 and _<m.shape[0]]
        cz = [_ for _ in cz if _>=0 and _<m.shape[0]]
        #print bead_i, cx[0], cx[-1], cy[0], cy[-1], cz[0], cz[-1]
        for i in cx:
            for j in cy:
                for k in cz:
                    m[i,j,k] = max(m[i,j,k] , tmp_sphere[i-cx[0], j-cy[0], k-cz[0]])
        del tmp_sphere, cx, cy, cz, t1, tmp
    IF.put_mrc(m, file_name)
    del m, rad, cords, rad_siz

def pool_function(farthest_point, c_com_i, inds, radii, scaling_factor, file_name_i):
    for struct_i in range(len(c_com_i)):
        file_name = os.path.join(file_name_i, ("chr12_" + str(inds[struct_i]) + ".mrc"))
        if os.path.isfile(file_name):   continue
        create_subtomogram(size=farthest_point, cords=c_com_i[struct_i], rad=radii, file_name=file_name, rescale=scaling_factor)
    return None

def main():
    return

if __name__=="__main__":
    main()
