#
#               Lab Exercise #1 for Lecture/Exercise 3.2
#                            Joe Felsenstein
#
#
# This is a lab exercise on contrasts, using the "pic" function of the
# "ape" phylogeny package in R.  I recommend that you *not* try to run
# this file, but instead keep it open in a window or an editor, and copy
# out the relevant parts and paste them into your R environment as needed.
#
# The data sets used will be 20-variable samples.  Each are
# measurements of the morphology of a fish, the fictional salmon-shark,
# whose form has evolved along a phylogeny.  The 20 variables come from
# the (x,y) coordinates of 10 landmarks on its two dimensional form:
# (x1,y1), (x2,y2), ... (x10, y10).  These forms have been superposed
# by a maximum likelihood method (Felsenstein and Bookstein, in prep.) and
# rescaled to remove size differences.  We will discuss that geometric
# morphometrics in a separate session.

# In this exercise, we use the word "coordinates" for the variables.
# If you are puzzled by the word coordinate, think of it as "variable".

# There is a tree, and the (imaginary) fishes evolve by correlated
# Brownian Motion, and the forms at the tips are sampled and superposed.
#
# There are 10 replicates of 100 fish in each sample.  The specimen names are
# A, AB, AC, B, BA, BB, BC, and so on to ZAB, but with name NA changed to
# NAA to avoid confusion with the R language symbol for missing information.

# There are five phylogenies and each has two replicate sets of fishes
# evolved on it.  Thus data sets fishes.1 and fishes.2 evolved in the first
# tree, fishes.3 and fishes.4 on the second, and so on. I have included
# those trees with the same numbering as the data sets.  So  intree.1
# is the tree for fishes.1, intree.2 is the tree for fishes.2, which is
# actually that same tree, intree.3 is the tree for fishes.3, which is
# the second tree, and so on.



# We start by getting the "ape" package
#
library(ape)
#

# If you don't have "ape" installed, this will not work.  In that case
# see the instructions at the Ape site:
#   http://ape-package.ird.fr/ape_installation.html
#
#
#
# Ape contains the functions read.tree and pic function that we will use.
# (pic  means "phylogenetically independent contrasts")
#
#
# (in these instructions I have been greatly helped by an excellent lab
# prepared by Brian O'Meara)
#
# which input file will we use:
#
infile = "fishes.1"
intree = "intree.1"
#
#
# Now read that tree into an  ape  tree
#
tr <- read.tree(intree)
#
#
# Now read the data into a frame 
#
datafile <- read.table(infile)
#
#
# how many species, how many coordinates
#
spp <- dim(datafile)[1]            # how many species
coords <- dim(datafile)[2] - 1     # how many columns not counting species names
#
# and also make a matrix with just the coordinates but without the species names
#
datacoords <- as.matrix(datafile[,2:(coords+1)])
#
#
# Now make a matrix with space for spp-1 contrasts per coordinate
#
contrasts <- matrix(0, spp-1, coords)
#
#
# Now things get obscure.  The "pic" function assumes that the rows of
# the data matrix are in the same order as the tips of the tree are, 
# going across the tree from left to right.
#
# The next statement reorders the rows of the data coordinates to make the rows
# of the data table (datacoords) be in that order:
#
datacoords <- datacoords[match(tr$tip.label,datafile[,1]),]
#
# NOTE:  this obscure step is VERY important, otherwise the mapping of tips
# to rows of the data matrix will be wrong, and you will get wierd,
# meaningless results.
#
#
#
# Now populate the columns by calling  pic  with tree  tr.  For each
# coordinate  pic (phylogenetically independent contrasts)  makes contrasts
# for one variable, which would be one column of our data matrix
#
for (i in 1:coords) { contrasts[,i] <- pic(datacoords[,i],tr) }
#
#
# Now we want to get covariances of characters, estimated using the
# contrasts and also a set estimated using the raw species values
#
contrastcovs <- cov(contrasts)
#
tipcovs <- cov(datacoords)
#
#
# It might be interesting to get the first principal component of each of
# these covariance matrices.
# In the R functions "eigen" and "svd" the first column of the matrix of
# eigenvectors are the weights of the first principal component
# (the second principal component, PC2, may be worth looking at too
#
contrastsPC1 <- eigen(contrastcovs)$vectors[,1]
#
tipsPC1 <- eigen(tipcovs)$vectors[,1]
#
#
# You can think of doing this on all the data sets  fishes.1, fishes.2, ...
# each with the corresponding tree (intree.1, intree.2, intree.3, ...
#
# If instead of the covariances we want the correlation coefficients
# we can use the R function "cov2cor":
#
contrastcorrs <- cov2cor(contrastcovs)
tipcorrs <- cov2cor(tipcovs)
#
# 
# If you do any plotting, it should be of tip values for one
# coordinate against tip values for another (i.e. one variable against
# another).  Or contrasts for one coordinate against contrasts for
# another.  Do NOT plot covariances against each other -- that way
# lies madness.
# 
# 
# Here is a question we can ask for simulation replicates but not for
# real data:  So, what is the truth?
# Answer:   Each coordinate underwent independent Brownian
# motion along the tree, except that coordinates 1, 3, 19 (variables V2, V4,
# and V20) also changed in a perfectly correlated Brownian motion that was
# somewhat larger, so that the nose evolved forward and backward,
# and in a second, independent set of correlated changes the dorsal fin
# went up and down (coordinate 7, variable V8) and the tip of the tail
# (coordinates 13 and 14 which are variables V14 and V15) changed in another
# Brownian motion, so that the top of the tail changed out and in along a
# line at a 45 degree angle..  There were also size changes which we hope the size
# correction has reduced in influence.
# 
# 
#    (For morphometrics data only)  Plotting the forms and PCs
#
# I have also included two R functions to plot the mean form of a dataset
# as a gray outline, and on it place red arrows showing the direction and
# magnitude of change for the two coordinates of each point in a
# principal component.
#
# These plotting functions and some functions they use are in file
# called    plotforms.R   which you can load into your R environment
# with this command
#
source("plotmeans.R")
#
# The plotting functions to use
# are:
#     plotmeansandpc(datafile, pc)
# and
#     plotmeansandpcmed(datafile, pc)
#
# where the first argument (a) is a data set (say  datacoords) and the second
# is one of the eigenvectors of the principal components, such as 
#
pc <- contrastsPC1
#
#  or
#
pc <- tipsPC1
#
# The second function,  plotmeansandpcmed   has the principal components
# adjusted to be sparser (this will be described in the lecture on
# morphometrics).
#
