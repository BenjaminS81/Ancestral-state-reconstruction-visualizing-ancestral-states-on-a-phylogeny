### Anole body size ancestral character state estimation simulation study.

### The script simulates estimates for ancestral states for a a continuous character, 
### in this case overall body size, on a given phylogeny of Anolis lizards from the Caribbean.

### To estimate the states for a continuously valued character at ancestral nodes, we need find
### the states that have the maximum probability of having arisen under our model. 
### These will be our maximum likelihood estimates. Ancestral character esimation is implemented
### in a variety of different R functions. The most commonly used is the ace function from Paradis.
### We start with the phytools function fastAnc. For these exercises we will use the phytools package, 
### as well as dependent packages such as ape. To follow all parts of the exercise you can 
### install the latest version of phytools - not from CRAN, but from the phytools page if needed.

###########################
### Start of simulation ###
###########################

## If phytools is loaded detach("package:phytools",unload=TRUE)
## Install dependencies for library ape if necessary. Check if R version is compliant with respective dependencies.
## It may require multiple R versions and switches in the global env. to install all packages correctly.

install.packages("http://www.phytools.org/nonstatic/phytools_0.4-98.tar.gz", type="source",repos=NULL)

## Installing package into 'C:/Users/Documents/R/win-library/3.2' (as 'lib' is unspecified)

packageVersion("phytools")

## Should return or higher: [1] '0.4.98'

## load library

library(phytools)

## Loading required package: ape
## Loading required package: maps

## Read tree from file.
## Note. The tree needs to be ultrametric with patristic distancace being equilibrated

anole.tree<-read.tree("anole.tre") ## or
anole.tree<-read.tree("http://www.phytools.org/eqg2015/data/anole.tre")

## Plot initial tree

plotTree(anole.tree,type="fan",ftype="i")

## setEnv=TRUE for this type is experimental. may have bugs with other types

## Read meta data

svl<-read.csv("svl.csv",row.names=1) ## or
svl<-read.csv("http://www.phytools.org/eqg2015/data/svl.csv", row.names=1)

## change this into a vector

svl<-as.matrix(svl)[,1]
svl

## Estimate ancestral states. Compute variances & 95% confidence intervals for each node. Call fit

fit<-fastAnc(anole.tree,svl,vars=TRUE,CI=TRUE)
fit

## Return 95% CI of 1, and return range of values within matrix

fit$CI[1,]
range(svl)

## Phytools has several different methods for visualizing reconstructed ancestral states 
## for a continuous trait on a given tree. One is a color gradient projection as below.

## Projection of the reconstruction onto the edges of the tree

obj<-contMap(anole.tree,svl,plot=FALSE)
plot(obj,type="fan",legend=0.7*max(nodeHeights(anole.tree)),fsize=c(0.7,0.9))

## A second projection of the tree into a space defined by time on the horizontal axis, and phenotype on the vertical dimension.

phenogram(anole.tree,svl,fsize=0.6,spread.costs=c(1,0))

## Properties of ancestral states of ancestral state reconstruction
## Next, we can explore some of the properties of ancestral state reconstruction of continuous traits in general.
## To do this, we will use some simulation functions in the phytools package.
## First, let's simulate some data:

## simulate a tree & some data

tree<-pbtree(n=26,scale=1,tip.label=LETTERS)

## simulate with ancestral states

x<-fastBM(tree,internal=TRUE)

## ancestral states

a<-x[as.character(1:tree$Nnode+Ntip(tree))]

## tip data

x<-x[tree$tip.label]

## Now, let's estimate ancestral states using fastAnc which uses the re-rooting algorithm discussed in class:

fit<-fastAnc(tree,x,CI=TRUE)
fit

## We can easily compare these estimates to the (normally unknown) generatied values for simulated data:

plot(a,fit$ace,xlab="true states",ylab="estimated states")

lines(range(c(x,a)),range(c(x,a)),lty="dashed",col="red") ## 1:1 line

## One of the most common critiques of ancestral state estimation is that the uncertainty around ancestral states can be large; and, 
## furthermore, that this uncertainty is often ignored. Let's address this by obtaining 95% CIs on ancestral values, and then let's 
## show that the 95% CIs include the generated values around 95% of the time:

## first, let's go back to our previous dataset accordingly

print(fit)

mean(((a>=fit$CI95[,1]) + (a<=fit$CI95[,2]))==2)

## One small tree doesn't tell us much, so let's repeat for a sample of trees & reconstructions:
## custom function that conducts a simulation, estimates ancestral states, & returns the fraction on 95% CI

foo<-function(){
    tree<-pbtree(n=100)
    x<-fastBM(tree,internal=TRUE)
    fit<-fastAnc(tree,x[1:length(tree$tip.label)],CI=TRUE)
    mean(((x[1:tree$Nnode+length(tree$tip.label)]>=fit$CI95[,1]) +
        (x[1:tree$Nnode+length(tree$tip.label)]<=fit$CI95[,2]))==2)
}

## conduct 100 simulations

pp<-replicate(100,foo())    
mean(pp)

## ## [1] 0.949899
## This shows us that although CIs can be large, when the model is correct they will include the generating value (1-α)% of the time.

## Ancestral state estimates when some nodes are known

## We can theoretically fix the state for any internal node during MLE of ancestral states. In practice, we would do this by attaching 
## a terminal edge of zero length to the internal node we wanted to fix, and then treat it like another tip value in our analyses.

## Another possibility, which also allows for the possibility that ancestral states are uncertain, is to use an informative prior probability 
## density on one or multiple states at internal nodes, and then estimate ancestral states using Bayesian MCMC.

## Let's try this using a really bad case for ancestral character estimation from contemporaneous tip data: Brownian evolution with a trend. 
## Note that although we theoretically could do so - we are not fitting a trended model here.

tree<-pbtree(n=100,scale=1)

## simulate data with a trend

x<-fastBM(tree,internal=TRUE,mu=3)
phenogram(tree,x,ftype="off",spread.labels=FALSE)
a<-x[as.character(1:tree$Nnode+Ntip(tree))]
x<-x[tree$tip.label]

## let's see how bad we do if we ignore the trend

plot(a,fastAnc(tree,x),xlab="true values",ylab="estimated states under BM")
lines(range(c(x,a)),range(c(x,a)),lty="dashed",col="red")
title("estimated without prior information")

## incorporate prior knowledge

pm<-setNames(c(1000,rep(0,tree$Nnode)), c("sig2",1:tree$Nnode+length(tree$tip.label)))

## the root & two randomly chosen nodes

nn<-as.character(c(length(tree$tip.label)+1,sample(2:tree$Nnode+length(tree$tip.label),2)))
pm[nn]<-a[as.character(nn)]

## prior variance

pv<-setNames(c(1000^2,rep(1000,length(pm)-1)),names(pm))
pv[as.character(nn)]<-1e-100

## run MCMC

mcmc<-anc.Bayes(tree,x,ngen=100000,
    control=list(pr.mean=pm,pr.var=pv,
    a=pm[as.character(length(tree$tip.label)+1)],
    y=pm[as.character(2:tree$Nnode+length(tree$tip.label))]))

## Control parameters (set by user or default):

## List of 7
##  $ sig2   : num 1.1
##  $ a      : Named num 0
##   ..- attr(*, "names")= chr "101"
##  $ y      : Named num [1:98] 0 0 0 0 0 ...
##   ..- attr(*, "names")= chr [1:98] "102" "103" "104" "105" ...
##  $ pr.mean: Named num [1:100] 1000 0 0 0 0 ...
##   ..- attr(*, "names")= chr [1:100] "sig2" "101" "102" "103" ...
##  $ pr.var : Named num [1:100] 1e+06 1e-100 1e+03 1e+03 1e+03 ...
##   ..- attr(*, "names")= chr [1:100] "sig2" "101" "102" "103" ...
##  $ prop   : num [1:100] 0.011 0.011 0.011 0.011 0.011 ...
##  $ sample : num 100
## Starting MCMC...
## Done MCMC.

nc.est<-colMeans(mcmc[201:1001,
    as.character(1:tree$Nnode+length(tree$tip.label))])
plot(a[names(anc.est)],anc.est,xlab="true values",
    ylab="estimated states using informative prior")
lines(range(c(x,a)),range(c(x,a)),lty="dashed",col="red")
title("estimated using informative prior")

## Discrete characters
## This part of the tutorial focuses on the estimation of ancestral character states for discretely valued traits using a continuous-time Markov chain model 
## commonly known as the Mk model.

## phytools" has a function for ancestral character estimation under a couple of different models that uses the re-rooting method of Yang (1995). 
## For other models of evolutionary change not covered by this function, users should try Rich Fitzjohn's diversitree package, or the ace function in ape.

## As a starting point pull the data from the Anolis ecomorph tree (packaged with phytools) to use in these analyses:

data(anoletree)

## this is just to pull out the tip states from the
## data object - normally we would read this from file

x<-getStates(anoletree,"tips")
tree<-anoletree
rm(anoletree)
tree

## Returns a phylogenetic tree with 82 tips and 81 internal nodes.
## Tip labels: Anolis_ahli, Anolis_allogus, Anolis_rubribarbus, Anolis_imias, Anolis_sagrei, Anolis_bremeri, ...
## Tree is rooted and includes branch lengths.
## Call x to depict tip labels

x

## In visual format, these are the data that we have:

plotTree(tree,type="fan",fsize=0.8,ftype="i")

## setEnv=TRUE for this type is experimental. please be patient with bugs

cols<-setNames(palette()[1:length(unique(x))],sort(unique(x)))
tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
    y=-max(nodeHeights(tree)),fsize=0.8)


## estimate ancestral states under a ER model
## fit a single-rate model & reconstruct ancestral states at internal nodes in the tree.

fitER<-ace(x,tree,model="ER",type="discrete")
fitER

## Example output from fitER:

## 
##     Ancestral Character Estimation
## 
## Call: ace(x = x, phy = tree, type = "discrete", model = "ER")
## 
##     Log-likelihood: -78.04604 
## 
## Rate index matrix:
##    CG GB TC TG Tr Tw
## CG  .  1  1  1  1  1
## GB  1  .  1  1  1  1
## TC  1  1  .  1  1  1
## TG  1  1  1  .  1  1
## Tr  1  1  1  1  .  1
## Tw  1  1  1  1  1  .
## 
## Parameter estimates:
##  rate index estimate std-err
##           1   0.0231   0.004
## 
## Scaled likelihoods at the root (type '...$lik.anc' to get them for all nodes):
##          CG          GB          TC          TG          Tr          Tw 
## 0.018197565 0.202238628 0.042841575 0.428609607 0.004383532 0.303729093

round(fitER$lik.anc,3)

## The element lik.anc gives us the marginal ancestral states, also known as the 'empirical Bayesian posterior probabilities.'
## It is fairly straightforward to overlay these posterior probabilities on the tree:

plotTree(tree,type="fan",fsize=0.8,ftype="i")

## setEnv=TRUE for this type is experimental. please be patient with bugs

nodelabels(node=1:tree$Nnode+Ntip(tree),
    pie=fitER$lik.anc,piecol=cols,cex=0.5)
tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
    y=-max(nodeHeights(tree)),fsize=0.8)

## An alternative approach to the one outline above is to use an MCMC approach to sample character histories from their posterior probability distribution. This is called stochastic character mapping (Huelsenbeck et al. 2003). The model is the same but in this case we get a sample of unambiguous histories for our discrete character's evolution on the tree - rather than a probability distribution for the character at nodes.
## For instance, given the data simulated above - we can generate the stochastic character map as follows:

## simulate single stochastic character map using empirical Bayes method

mtree<-make.simmap(tree,x,model="ER")

## Example output from simmap:
## make.simmap is sampling character histories conditioned on the transition matrix
## Q =
##             CG          GB          TC          TG          Tr          Tw
## CG -0.11570723  0.02314145  0.02314145  0.02314145  0.02314145  0.02314145
## GB  0.02314145 -0.11570723  0.02314145  0.02314145  0.02314145  0.02314145
## TC  0.02314145  0.02314145 -0.11570723  0.02314145  0.02314145  0.02314145
## TG  0.02314145  0.02314145  0.02314145 -0.11570723  0.02314145  0.02314145
## Tr  0.02314145  0.02314145  0.02314145  0.02314145 -0.11570723  0.02314145
## Tw  0.02314145  0.02314145  0.02314145  0.02314145  0.02314145 -0.11570723
## (estimated using likelihood);
## and (mean) root node prior probabilities
## pi =
##        CG        GB        TC        TG        Tr        Tw 
## 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667
## Done.

mtree

## Example output:

## Phylogenetic tree with 82 tips and 81 internal nodes.
## 
## Tip labels:
##  Anolis_ahli, Anolis_allogus, Anolis_rubribarbus, Anolis_imias, Anolis_sagrei, Anolis_bremeri, ...
## 
## The tree includes a mapped, 6-state discrete character with states:
##  CG, GB, TC, TG, Tr, Tw
## 
## Rooted; includes branch lengths.

plot(mtree,cols,type="fan",fsize=0.8,ftype="i")

## setEnv=TRUE for this type is experimental. please be patient with bugs

add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
    y=-max(nodeHeights(tree)),fsize=0.8)

## A single stochastic character map does not mean a whole lot in isolation - we need to look at the whole distribution 
## from a sample of stochastic maps. This can be a bit overwhelming. For instance, the following code generates 100 stochastic 
## character maps from our dataset and plots them in a grid:

mtrees<-make.simmap(tree,x,model="ER",nsim=100)

## Example output:
## make.simmap is sampling character histories conditioned on the transition matrix
## Q =
##             CG          GB          TC          TG          Tr          Tw
## CG -0.11570723  0.02314145  0.02314145  0.02314145  0.02314145  0.02314145
## GB  0.02314145 -0.11570723  0.02314145  0.02314145  0.02314145  0.02314145
## TC  0.02314145  0.02314145 -0.11570723  0.02314145  0.02314145  0.02314145
## TG  0.02314145  0.02314145  0.02314145 -0.11570723  0.02314145  0.02314145
## Tr  0.02314145  0.02314145  0.02314145  0.02314145 -0.11570723  0.02314145
## Tw  0.02314145  0.02314145  0.02314145  0.02314145  0.02314145 -0.11570723
## (estimated using likelihood);
## and (mean) root node prior probabilities
## pi =
##        CG        GB        TC        TG        Tr        Tw 
## 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667
## Done.

mtrees

## 100 phylogenetic trees with mapped discrete characters

par(mfrow=c(10,10))
null<-sapply(mtrees,plot,colors=cols,lwd=1,ftype="off")

## It's possible to summarize a set of stochastic maps in a much more meaningful way. For instance, we can estimate the 
## number of changes of each type, the proportion of time spent in each state, and the posterior probabilities that each 
## internal node is in each state, under our model. For example:

pd<-summary(mtrees,plot=FALSE)

pd

## 100 trees with a mapped discrete character with states as example output:
##  CG, GB, TC, TG, Tr, Tw 
## 
## trees have 24.07 changes between states on average
## 
## changes are of the following types:
##      CG,GB CG,TC CG,TG CG,Tr CG,Tw GB,CG GB,TC GB,TG GB,Tr GB,Tw TC,CG
## x->y  0.27  0.25  0.22  0.09  0.23  0.56  0.97  1.26  0.41  1.52  1.41
##      TC,GB TC,TG TC,Tr TC,Tw TG,CG TG,GB TG,TC TG,Tr TG,Tw Tr,CG Tr,GB
## x->y  0.64  0.32   0.8  1.13  0.97  3.17  1.91  0.93  1.81   0.2  0.17
##      Tr,TC Tr,TG Tr,Tw Tw,CG Tw,GB Tw,TC Tw,TG Tw,Tr
## x->y  0.15  0.14  0.13  0.56  1.31  1.43   0.7  0.41
## 
## mean total time spent in each state is:
##               CG        GB         TC         TG          Tr         Tw
## raw  12.90191857 47.001692 33.6593396 69.1939136 12.53922339 30.3708734
## prop  0.06273209  0.228533  0.1636594  0.3364367  0.06096858  0.1476702
##        total
## raw  205.667
## prop   1.000

plot(pd,fsize=0.6,ftype="i")

## now let's plot a random map, and overlay the posterior probabilities

plot(mtrees[[1]],cols,type="fan",fsize=0.8,ftype="i")

## setEnv=TRUE for this type is experimental. please be patient with bugs

nodelabels(pie=pd$ace,piecol=cols,cex=0.5)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
    y=-max(nodeHeights(tree)),fsize=0.8)

## Finally, since we obtained these inferences under exactly the same model, let's compare the posterior probabilities from stochastic mapping
## with our marginal ancestral states. In the former case, our probabilities were obtained by sampling from the joint (rather than marginal) 
## probability distribution for the ancestral states.

plot(fitER$lik.anc,pd$ace,xlab="marginal ancestral states",
    ylab="posterior probabilities from stochastic mapping")
lines(c(0,1),c(0,1),lty="dashed",col="red",lwd=2)

## This tells us that although joint & marginal reconstruction are not the same, the marginal probabilities from joint stochastic sampling 
## and the marginal ancestral states are quite highly correlated - at least in this case study.

## Other visualization methods for ancestral states:

## In this final part of the tutorial I'm going to quickly overview a range of plotting methods for phylogenies & comparative data 
## that are implemented in the phytools package.

## Continuous character methods
## The first plotting method I'll illustrate is called a 'traitgram' which is a projection of the tree into a space defined by trait 
## values & time since the root. So, for example:

## First simulate a tree

tree<-pbtree(n=26,tip.label=LETTERS)
plotTree(tree)

## simulate data under Brownian motion

x<-fastBM(tree)
x

## Plot the BM phenogram

phenogram(tree,x,spread.labels=TRUE,spread.cost=c(1,0))

## We can also plot a traitgram with the uncertainty about ancestral character states visualized using transparent probability density. 
## This is implemented in the phytools function fancyTree. So, for instance:

## plot traitgram with 95% CI

fancyTree(tree,type="phenogram95",x=x,spread.cost=c(1,0))

## Next, we can explore the continuous character mapping function we've seen already in phytools called contMap. Here is a quick demo using the same data:

## plot contMap

obj<-contMap(tree,x)

## The function contMap returns an object of class "contMap" which we can then more easily replot using, for instance, different parameters:

## plot leftward

plot(obj,direction="leftwards")

## For large trees, we might want to use a circular or “fan” tree:

tree<-pbtree(n=200,scale=1)
x<-fastBM(tree)
obj<-contMap(tree,x,plot=FALSE)
plot(obj,type="fan",outline=FALSE,ftype="off")

## Note that the intense aliasing is just a byproduct of plotting to .png format. To get a high quality image we can plot to .pdf or other formats. For instance:

pdf(file="plot-contMap.pdf",width=10,height=10)
plot(obj,type="fan",outline=FALSE)
dev.off()

## Finally, a relatively simple new plotting method in phytools is the function plotTree.wBars. That function pretty much does what it sounds like it does:

plotTree.wBars(anole.tree,exp(svl),type="fan",scale=0.002,
    fsize=0.7,tip.labels=TRUE)

## It is not too difficult to combine this with a contMap plot. For example:

obj<-contMap(anole.tree,exp(svl),plot=FALSE)
plotTree.wBars(obj$tree,exp(svl),method="plotSimmap",
    tip.labels=TRUE,fsize=0.7,colors=obj$cols,type="fan",scale=0.002)
add.color.bar(1.0,obj$cols,title="trait value",lims=obj$lims,prompt=FALSE,
    x=0.9*par()$usr[1],y=0.9*par()$usr[3])

## We can also do two dimensional visualizations of phylogenies in morphospace. The is called a 'phylomorphospace'. E.g.:

tree<-pbtree(n=26,tip.label=LETTERS)
X<-fastBM(tree,nsim=3) ## simulate 3 characters
colnames(X)<-paste("trait",1:3)
phylomorphospace(tree,X[,c(1,2)],xlab="trait 1",ylab="trait 2")

## With phytools we can do static or animated version of the same. I'll just do a static version here, but animated is the default:

phylomorphospace3d(tree,X,method="static")

## Finally, it's possible to combine the continuous character mapping and 2D phylomorphospaces using a type of phylogenetic scatterplot. Here's a demo:

fancyTree(tree,type="scattergram",X=X)

## Discrete character methods
## Stochastic character maps are easy to plot using phytools. For instance:

data(anoletree)

plotSimmap(anoletree,fsize=0.6,ftype="i",ylim=c(-1,Ntip(anoletree)))

## no colors provided. using the following legend:
##        CG        GB        TC        TG        Tr        Tw 
##   "black"     "red"  "green3"    "blue"    "cyan" "magenta"
## we can use any color scheme, but this is the default

cols<-setNames(palette()[1:length(unique(getStates(anoletree,
"tips")))],sort(unique(getStates(anoletree,"tips"))))
add.simmap.legend(colors=cols,x=0,y=-1,vertical=FALSE,prompt=FALSE)

## mark the changes on the tree

xy<-markChanges(anoletree,plot=FALSE)
points(xy,pch=19)

## For binary characters, we can also map the posterior density of stochastic maps onto the nodes & branches of a 
## tree using the phytools function densityMap. Here's a demo of this using simulated data.

## first let's simulate data

Q<-matrix(c(-1,1,1,-1),2,2)
rownames(Q)<-colnames(Q)<-c(0,1)
tree<-sim.history(tree,Q)
x<-tree$states
x

##   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R 
## "1" "1" "1" "0" "0" "1" "1" "1" "1" "1" "0" "1" "0" "1" "1" "0" "0" "0" 
##   S   T   U   V   W   X   Y   Z 
## "1" "1" "1" "1" "0" "0" "0" "0"

## stochastic maps via simmap

trees<-make.simmap(tree,x,nsim=100)

## Example output:

## make.simmap is sampling character histories conditioned on the transition matrix
## Q =
##            0          1
## 0 -0.9914265  0.9914265
## 1  0.9914265 -0.9914265
## (estimated using likelihood);
## and (mean) root node prior probabilities
## pi =
##   0   1 
## 0.5 0.5
## Done.
## make densityMap

obj<-densityMap(trees,lwd=4,outline=TRUE)

## Combining discrete & continuous characters
## We can, for instance, overlay a discrete character history (or stochastic map) onto a traitgram or phylomorphospace.
## Simulate data with a high rate on some branches as follows:

x<-sim.rates(tree,setNames(c(1,10),c(0,1)))

## Plot traitgram

phenogram(tree,x,colors=setNames(c("blue","red"),
    c(0,1)),spread.labels=TRUE,spread.cost=c(1,0))

## Plot phylomorphospaces

X<-cbind(sim.rates(tree,setNames(c(1,10),c(0,1))),
    sim.rates(tree,setNames(c(1,10),c(0,1))))
phylomorphospace(tree,X,colors=setNames(c("blue","red"),c(0,1)))

## Overlay a posterior density from stochastic mapping

phylomorphospace(obj$tree,X,colors=obj$cols,
    lwd=3,xlab="x",ylab="y")
add.color.bar(4,obj$cols,title="PP(state=1)",
    prompt=FALSE,x=0.9*par()$usr[1],y=0.9*par()$usr[3])

#####################################################################################
## End of simulation and case study. If you have any questions, please let me know ##    
#####################################################################################


