#Here is a script to help you to obtain the estimates of evolutionary rates for discrte traits and simulate traits under a given model.
library(geiger)
library(phytools)
#Since I cannot at this moment access the tree and data we have, I will use simulated data, so this next part of the code is simply to simulate a 
#tree and data to work with, you won't need this for your analses.

tree<-sim.bdtree(b=1, d=0, stop="taxa", n=100, extinct=FALSE)

#Now simulate some data

#For discrete data, we need a Q matrix, which simply describes the transition matrix between states. You will see these transition
#matrices lower down when comparing different models of trait evolution for discrete characters, they are square matrices where by definition
#the sum of row values are equal to zero. For two states of a trait (0,1) the matrix will be 2x2, for 3 states it will be 3x3, and so on.

Q<-matrix(c(-0.40,0.30,0.1,0.5,-0.70,0.2,0.15,0.15,-0.3), nrow=3, ncol=3, byrow=TRUE)
data<-as.data.frame(sim.char(tree, par=Q, nsim=2, model="discrete"))
colnames(data)<-c("Trait1", "Trait2")

#You see now we have a data set with 2 discrete traits with 3 states.

##############################From here this code you will need#################################################

#First off we need to estimate the model of trait evoltuion:
#This fits an "equal rates" model, where all transitions have equal value
trait1<-data$Trait1
names(trait1)<-row.names(data)

get_traitQ <- function(trait, tree) {
  trait.er<-rerootingMethod(tree, trait, model="ER")
  trait.sym<-rerootingMethod(tree, trait, model="SYM")
  if(pchisq(-2*(trait.er$loglik-trait.sym$loglik), df=3-1, lower.tail=FALSE) > 0.05){
    Q <- trait.er$Q
    attr(Q, which = "model") <- "ER"
  }else{
    Q <- trait.sym$Q
    attr(Q, which = "model") <- "SYM"
  }
  return(Q)
}

get_traitSim <- function(Q, tree, nsim = 500) {
  sims <- as.data.frame(sim.char(tree, par = Q, nsim = nsim, model="discrete"))
  return(sims)
}

get_simGKtau <- function(pair, data, tree, nsim = 500) {
  
sims1 <- get_traitSim(Q = get_traitQ(traits[,1], tree), tree, nsim = nsim)

sims2 <- get_traitSim(Q = get_traitQ(traits[,2], tree), tree, nsim = nsim)


tau <- mapply(FUN = function(x, y){
  tau <- GKtau(x,y)
  tau[, c("tauxy", "tauyx")]}, x = sims1, y = sims2, SIMPLIFY = F) %>% do.call(what = rbind)

return(tau)
}



colMeans(tau)

GKtau(x = trait1, y = trait2)

trait1.er<-rerootingMethod(tree, trait1, model="ER")
#This fits a "symetric model" where forward and backward transitions between given states take on the same value, but transitions between
#different states can differ (so e.g. 1 to 2, takes on the same value as 2 to 1, but 2 to 3 and 3 to 2 can take on a different value than 1 to 2)
trait1.sym<-rerootingMethod(tree, trait1, model="SYM")

#for traits with only two states use "ARD" instead of "SYM", but I would not recommend using ARD for traits with many states because it becomes
#difficult to fit and you can get overestimaes of goodness-of-fit (likelihood) for this model

#Compare the fit of the different models:

trait1.er$loglik
trait1.sym$loglik

# Compare log likelihood between the models 
# Log likelihood ratio test = -2 * Log likelihood of the simpler model minus log 
# lik of the more complex model, df = number of parameters in the more complex model 
# minus the number of parameters in the simpler model. The more complex model always 
# has a likelihood that is the same or higher than the simpler model, if not there 
# is a problem with model fit. 

pchisq(-2*(trait1.er$loglik-trait1.sym$loglik), df=3-1, lower.tail=FALSE)

#In this case there is not a significant difference between the ER and SYM models and 
# thus we chose the ER model (because of parsimony) We can now use this model, 
# with the estimated Q matrix to simulate a bunch of Traits1 under a similar model
# of evolution as trait1.

Qt1<-trait1.er$Q
Trait1.syms<-as.data.frame(sim.char(tree, par= Qt1, nsim=500, model="discrete"))

#We can do the same for Trait2
trait2<-data$Trait2
names(trait2)<-row.names(data)
trait2.er<-rerootingMethod(tree, trait2, model="ER")
trait2.sym<-rerootingMethod(tree, trait2, model="SYM")
#Compare log likelihood between the models 
pchisq(-2*(trait2.er$loglik-trait2.sym$loglik), df=3-1, lower.tail=FALSE)

Qt2<-trait2.er$Q
Trait2.syms<-as.data.frame(sim.char(tree, par= Qt2, nsim=500, model="discrete"))

#Now you simply recreate a data frame with 500 (or whatever number is justified) 
#simulated variables 1 and 500 simulated variables 2, which since they are simulated 
#independently should not be correlated, with which you can estimate a distribution 
#of expected GoodmanKruskal statistics against which to compare actual values. 
#Hope this is clear! 
