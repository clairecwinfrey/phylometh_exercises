##Claire Winfrey
#PhyloMeth
#April 3,2020
#Diversification


##############################################################################################################
#                                          TREE ONLY MODELS 
##############################################################################################################

#install.packages(c("ape", "TreeSim", "geiger", "diversitree", "devtools"))
library(ape)
library(TreeSim)
library(geiger)
library(diversitree)
#devtools::install_github("thej022214/hisse")
library(hisse)

#SIMULATE A TREE WITH NO EXTINCTION with 300 taxa, only speciation, no extinction
my.tree <- TreeSim::sim.bd.taxa(n=300, numbsim=1, lambda=0.1, mu=0)[[1]]

#Now, let's plot the tree using a lineage through time plot:

quartz()
ape::ltt.plot(my.tree)
#Shows that species are increasing exponentially over time.

#Put it on a log scale:
ape::ltt.plot(my.tree, log ="y")

#LOOKING AT MULTIPLE TREES
#sim.bd.taxa function takes a fixed number of taxa and simulates birth and death trees (constant process). 
#complete option lets you decide if you want to have tree with or without extinct lineages.
#lambda is speciation rate, mu is extinction rate
yule.trees <- TreeSim::sim.bd.taxa(n=300, numbsim = 10, lambda = 0.1, mu = 0, complete = FALSE)
#to plot multiple trees:
ape::mltt.plot(yule.trees, log= "y")

#TREES WITH BIRTH AND DEATH (just change mu parameter to non-zero)

bd.trees <- TreeSim::sim.bd.taxa(n=300, numbsim = 10, lambda = 1, mu =0.9, complete = FALSE)
quartz()
ape::mltt.plot(bd.trees, log= "y", legend =FALSE)

#COMPARE TREES WITH ONLY BIRTH VERSUS BIRTH AND DEATH
depth.range <- range(unlist(lapply(yule.trees,ape::branching.times)), unlist(lapply(bd.trees,ape::branching.times)))
# 1) range returns a vector containing the minimum and maximum of all the given arguments 2) Unlist turns list into a vector. 
#3) lapply returns a list of the same length as x (lenght of yule trees is 10 trees). #4) Branching times computes distance 
#from each node to the tips of a tree 

#So this code returns two numbers, which are the branching times for the birth only and the birth-death tree, respectively.

max.depth <- sum(abs(depth.range)) #ape rescales depths
plot(x=c(0, -1*max.depth), y=c(1, ape::Ntip(yule.trees[[1]])), log="y", type="n", bty="n", xlab="Time", ylab="N")
colors=c(rgb(1,0,0,0.5), rgb(0, 0, 0, 0.5))
list.of.both <- list(bd.trees, yule.trees)

for (i in sequence(2)) {
  tree.list <- list.of.both[[i]]
  for (j in sequence(length(tree.list))) {
    ape::ltt.lines(tree.list[[j]], col=colors[[i]])
  }
}
legend("topleft", legend=c("Birth Death", "Yule"), fill=colors)

#zoom in on final part of plot:
depth.range <- range(unlist(lapply(yule.trees,ape::branching.times)), unlist(lapply(bd.trees,ape::branching.times)))
max.depth <- sum(abs(depth.range)) #ape rescales depths
plot(x=c(0, -5), y=c(200, ape::Ntip(yule.trees[[1]])), log="y", type="n", bty="n", xlab="Time", ylab="N")
colors=c(rgb(1,0,0,0.5), rgb(0, 0, 0, 0.5))
list.of.both <- list(bd.trees, yule.trees)
for (i in sequence(2)) {
  tree.list <- list.of.both[[i]]
  for (j in sequence(length(tree.list))) {
    ape::ltt.lines(tree.list[[j]], col=colors[[i]])
  }
}
legend("topleft", legend=c("Birth Death", "Yule"), fill=colors)

#EXPERIMENTING WITH OTHER DIVERSIFICATION PATTERNS

### 1. speciation rate much higher than extinction rate
my.trees1 <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda = 1, mu= 0.5, complete=FALSE)
quartz()
ape::mltt.plot(my.trees1, log="y", legend=FALSE)

#The curve starts shooting 

### 2. constant difference between lambda and mu, but different values
#lower values
my.trees2a <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda = .75, mu= 0.25, complete=FALSE)
ape::mltt.plot(my.trees2a, log="y", legend=FALSE)
?sim.bd.taxa()

#higher values
my.trees2b <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda = .75, mu= 0.25, complete=FALSE)
ape::mltt.plot(my.trees2b, log="y", legend=FALSE)

### 3. sum of lambda and mu is same, but different values!
#same values of lambda and mu
my.trees4a <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda = 0.75, mu= 0.75, complete=FALSE)
ape::mltt.plot(my.trees4a, log="y", legend=FALSE)

#less difference between lambda and mu
my.trees4b <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda = 1.25, mu= 0.25, complete=FALSE)
ape::mltt.plot(my.trees4b, log="y", legend=FALSE)

