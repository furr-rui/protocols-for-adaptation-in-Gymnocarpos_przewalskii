## sigma mean with uniform dist
## alpha mean with uniform dist
## theta mean with uniform dist 
#alpha.mean = seq(0.1,10,0.01)
set.seed()
alpha.mean = 10^seq(-2,1,0.01)

#sigma.mean = seq(0, 10,0.01)
sigma.mean = 10^seq(-2, 1,0.01)

#theta.mean = seq(0, 10, 0.01)
theta.mean = 10^seq(-2, 1,0.01)

num.genes.scale = seq(3,6,0.001)

n = 1e6
sigma.prior=sample(sigma.mean, size=n, rep=T)
alpha.prior=sample(alpha.mean, size=n, rep=T)
theta.prior=sample(theta.mean, size=n, rep=T)
num.genes=sample(num.genes.scale, size=n, rep=T)*1e4

write.table(data.frame(num.genes=num.genes, sigma=sigma.prior, alpha=alpha.prior, theta = theta.prior),file='~/Desktop/ZJU/Gymnocarpos/correlation/abc/prior/parameters_priors.txt',quote=F,row.names=F,sep="\t")