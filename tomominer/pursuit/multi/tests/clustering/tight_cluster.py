
'''

# code for testing performing tight clustering through the method of Tseng05

'''




'''
# R test code

rm(list = ls(all=TRUE)) 

n_norm = 50
n_unif = 1000
n_dim = 10

#Sigma = matrix(c(10,3,3,2),2,2)

library(MASS)
x_norm = mvrnorm(n=n_norm, mu=rep(0, n_dim), Sigma=diag(0.01, nrow=n_dim))
x_unif = (matrix(runif(n_unif * n_dim), n_unif, n_dim) - 0.5) * 40

x = rbind(x_norm, x_unif)
dim(x)

library(tightClust)
tclust1 = tight.clust(x, target=5, k.min=10, random.seed=12345)


tclust1$cluster[1:100]

xd = density(x, bw='SJ')
xd$y[1:100]


hc = hclust(dist(x))
plot(hc)

'''


