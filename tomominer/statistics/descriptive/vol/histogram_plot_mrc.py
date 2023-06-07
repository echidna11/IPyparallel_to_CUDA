#!/usr/bin/env python

'''
given a list of mrc files, plot the histogram for each
'''



if __name__ == '__main__':

    import matplotlib
    matplotlib.use('Qt4Agg')

    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt
       
    import sys
    import tomominer.io.file as IF
    for i, fn in enumerate(sys.argv):
        if i == 0:  continue

        v = IF.read_mrc_vol(fn)
        print 'mean', v.mean(), 'std', v.std()
        n, bins, patches = plt.hist(v.flatten(), 1000)
        plt.title(fn)
        plt.show()

