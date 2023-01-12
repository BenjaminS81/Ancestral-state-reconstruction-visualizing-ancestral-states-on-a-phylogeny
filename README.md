
## Ancestral character state estimate simulation study for continuous and discrete characters ##

In this highly guided R studio script the principles of infering ancestral character states employing different
mathematical models are simulated and discussed. The script mainly employs functions from the R packages Phytools and ape.

It can be in principle easily adopoted to any type of dataset and comparative quantiative biological question.

Files in this repo:

ACE.r : R code and simulation study
.tre : Tree file with Anole phylogeny (see below)
.csv : Metadata with body size data (see below)

A short description of the simulation:

Anole body size ancestral character state estimation simulation study.

The script simulates estimates for ancestral states for a a continuous character, 
in this case overall body size, on a given phylogeny of Anolis lizards from the Caribbean.

To estimate the states for a continuously valued character at ancestral nodes, we need find
the states that have the maximum probability of having arisen under the model. 
These will be the maximum likelihood estimates. Ancestral character esimation is implemented
in a variety of different R functions. The most commonly used is the ace function from Paradis.
We start with the phytools function fastAnc. For these exercises we will use the phytools package, 
as well as dependent packages such as ape. To follow all parts of the exercise you can 
install the latest version of phytools - not from CRAN, but from the phytools page if needed.
