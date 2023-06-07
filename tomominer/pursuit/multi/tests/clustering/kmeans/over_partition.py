'''
given data of one cluster, use kmeans with k=2 to see how the clustering behave
'''



import matplotlib
matplotlib.use('Qt4Agg')



import numpy as N

x = N.random.multivariate_normal([0,0], [[1,0], [0,1]], size=100)

import sklearn.cluster as SC
km = SC.KMeans(n_clusters=3)
l = km.fit_predict(x)


import matplotlib.pyplot as plt
plt.scatter(x[l==0,0], x[l==0,1], c='b')
plt.scatter(x[l==1,0], x[l==1,1], c='r')
plt.scatter(x[l==2,0], x[l==2,1], c='y')


plt.show()



#--------------------------
# two real clusters, kmeans k=3

import numpy as N
x0 = N.random.multivariate_normal([-5,0], [[1,0], [0,1]], size=100)
x1 = N.random.multivariate_normal([5,0], [[1,0], [0,1]], size=100)

x = N.vstack((x0, x1))

import sklearn.cluster as SC
km = SC.KMeans(n_clusters=3)
l = km.fit_predict(x)


import matplotlib.pyplot as plt
plt.scatter(x[l==0,0], x[l==0,1], c='b')
plt.scatter(x[l==1,0], x[l==1,1], c='r')
plt.scatter(x[l==2,0], x[l==2,1], c='y')


plt.show()



