rm(list=ls())


library(caper)
library(geiger)
library(PHYLOGR)

# brownian case, lambda = 1
nsps <- 300
tree <- sim.bdtree(n=nsps) 
x<- scale(rTraitCont(rescale(tree, "lambda", 1)))
y<- scale(rTraitCont(rescale(tree, "lambda", 1)))
x.pic<- pic(x, tree)
y.pic<- pic(y, tree)
cor.origin(x.pic,y.pic) # this will equal the cor.xy below

invC <-solve(vcv.phylo(tree))
mean.x <-colSums(invC%*%x)/sum(invC)
mean.y <-colSums(invC%*%y)/sum(invC)
vector.ones<-as.matrix(rep(1,nsps))
var.x <-t(x-vector.ones%*%mean.x) %*% invC%*%(x-vector.ones%*%mean.x)/(nsps-1)
var.y <-t(y-vector.ones%*%mean.x) %*% invC%*%(y-vector.ones%*%mean.x)/(nsps-1)
cor.xy <-(t(x-vector.ones%*%mean.x) %*% invC%*%(y-vector.ones%*%mean.x)/(nsps-1))/sqrt(var.x*var.y)
cor.xy

# brownian case, lambda <> 1
nsps <- 300
tree <- sim.bdtree(n=nsps) 
x<- scale(rTraitCont(rescale(tree, "lambda", 0.7)))
y<- scale(rTraitCont(rescale(tree, "lambda", 0.7)))

d<- data.frame(spp=tree$tip.label, x=x, y=y)
dat<- comparative.data(tree, d, 'spp', vcv=TRUE)
result.pgls<-pgls(y ~ x, dat,lambda="ML")

# use estimated lambda to rescale the covariance matrix
invC<-solve(vcv.phylo(rescale(tree, "lambda",result.pgls$param[2])))

mean.x<-colSums(invC%*%x)/sum(invC)
mean.y<-colSums(invC%*%y)/sum(invC)
vector.ones<-as.matrix(rep(1,nsps))
var.x<-t(x-vector.ones%*%mean.x)%*%invC%*%(x-vector.ones%*%mean.x)/(nsps-1)
var.y<-t(y-vector.ones%*%mean.x)%*%invC%*%(y-vector.ones%*%mean.x)/(nsps-1)
cor.xy<-(t(x-vector.ones%*%mean.x)%*%invC%*%(y-vector.ones%*%mean.x)/(nsps-1))/sqrt(var.x*var.y)
cor.xy

