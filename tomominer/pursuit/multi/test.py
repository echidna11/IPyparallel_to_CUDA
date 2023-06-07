import tomominer.io.file as IF
import numpy as N

def pool_function(farthest_point, c_com_i, inds, radii, scaling_factor):
    for i in range(len(inds)):
        tmp = N.zeros((4,4,4))
        file_name = str(inds[i]) + ".mrc"
        IF.put_mrc(tmp, file_name)

def main():
    return

if __name__=="__main__":
    main()
