# R exercise for Wed June 12th Quantitative Genetics workshop 2019
# Samantha Price 

# Section 1 covers basic Brownian motion on a phylogeny
# Section 2 covers changes in Brownian motion rate on a phylogeny 
## Part A changes over time - using PIC and model fitting 
## Part B changes over branches using model fitting

library(ape)
library(phytools)
library(geiger)
library(OUwie)

# For section 2 you will need the dataset of cetacean body mass and female length at sexual maturity in CetaceaLH.txt and the tree CetaceaPhy.nex

#==============SECTION 1 Brownian motion (BM) on a phylogeny============== 

# What is Brownian motion?  At each instant in time the state of a character can increase or decrease and the direction and magnitude of these changes are independent of current or past states (no trends etc.).  The total displacement of a character after time t is drawn from a normal distribtion with a mean of zero (i.e. net change in trait is 0) and a variance proportional to t. 

# Generate a basic 3 species tree using the ape read.tree function but you are giving it a newick string directly as text rather than a file

tree<-read.tree(text="((grizzlybear:1, polarbear:1):4, spectacledbear:5);")
plot(tree)
axisPhylo() # this plots an axis on the Phylogeny
nodelabels() # this plots the nodelabels, it tells you how ape labels the nodes

# Looking at this tree we see that grizzly bears and polar bears share about 4 million years of history during which they all follow the same BM trajectory after which grizzly and polar bears split with both lineages starting on a new BM trajectory - we can simulate this. We need the branch length information: the longer the branch the greater the expected variance - with this we determine the standard deviation of the normal distribution. The total displacement of a character after time v is drawn from a normal distribtion with a mean of zero (i.e. net change in trait is 0) and a standard deviation = sqrt(rate*branchlength). 

# Lets assume we have a root state of 1 and a sigma^2 of 0.5 

# Start with the first branch from the root node 4 to 5 and estimate the trait value at node 5 and use the branch lengths to estimate the displacement.

x<-rnorm(1, mean=0, sd=sqrt(0.5*tree$edge.length[1]))

# You obviously need to convert this displacement to a trajectory, so if the root state is 1 then you add the displacement to the root state to get the value at the node.

node5est<-1+x

# You then repeat this process across the tree calculating the displacement with sqrt(0.5*tree$edge.length[i]) and adding it to the ancestral node. Fortunately you don't need to be able to do this as there are many R functions that will do it for you, such as fastBM in phytools.

# We will now use simulations to explore the stochastic nature of the BM process and how differences in rate influence tip variance - feel free to use whatever methods you are familiar with but this exercise will use the phytools and the function phenogram to visualize the variance and the function fastBM to simulate traits with different rates of BM. 

# First simulate a phylogeny, my example will use 30 species birth rate of 1 and death rate of 0.02 stopping at 30 taxa using the rphylo function in ape

simtr<-rphylo(30, b=1, d=0.02, fossils=FALSE)

# Now simulate 4 traits with two different rate categories, feel free to pick your own, this example we will use BM rate = 1 and BM rate = 2. FastBM generates a matrix with nrows=number of tips and ncols=nsim

bm1sim<-fastBM(simtr, sig2=1, nsim=4)

bm2sim<-fastBM(simtr, sig2=2, nsim=4)


# Now visualize the trait evolution using phenograms which plots a projection of the phylogenetic tree in a space defined by phenotype on the y axis and time on the x axis.

par(mfrow=c(2,4))
for(i in 1:4) {phenogram(simtr, x= bm1sim[,i], col="grey")} # this short loop just saves us from manually running the phenogram function x4 on each of the different simulations
for(i in 1:4) {phenogram(simtr, x= bm2sim[,i], col="blue")}

# How does the variation within a rate category compare to between the two rate categories?
# Do you think you would be able to distinguish if a trait has evolved at with a rate of 1 or 2 on this phylogeny? Does your neighbour come to the same conclusion?

#==============SECTION 2 part A Changes in BM rate over time ============== 

cetphy<-read.nexus("CetaceaPhy.nex") # Cetacean phylogeny 

cetdat<-read.table("CetaceaLH.txt", header=T) #Cetacean life history data

tmp<-treedata(cetphy, cetdat, sort=T) # useful function in geiger that creates an object with the overlapping taxa in the tree and the dataset

tr<-tmp[[1]]
dat<-tmp[[2]]

# We will first use the node height test on independent contrasts to determine if the rate of Female length at sexual maturity is speeding up or slowing down over time 

# First estimate the log of mass which is in the first column,

logflsm<-log10(dat[,2])

# Then calculate the node height, you could do this in several steps, using the pic function in ape and then fitting a linear model on the absolute value of the contrasts calculated and height above the root of the node at which they were being compared but there is a function in geiger: nh.test that does it all for you

logflsmNH<-nh.test(tr, logflsm, regression.type="lm", show.plot=TRUE) # the other option is a robust regression, which can be very useful see Slater & Pennell 2013 syst. biol.

# Why does Brownian motion specifically predict that the absolute value of the standardize contrast should not be related to age? So what does it mean if the absolute contrasts of Female length at sexual maturity are higher closer to the root? 

# We will now use model-fitting to determine if Female length at sexual maturity is best fit by a Brownian motion model, a model that allows the rate to speed up or slow down over time or a non-phylogenetic white noise model that assumes a single normal distribution of the data with no phylogenetic covariance amongst specie. We will use the fitContinuous function in geiger for these models and AIC to identify the best-fitting model (problems with AIC are likely to be discussed later). When using the Akaike Information Criterion (corrected for small sample size AICc) the rule of thumb is a difference between the best fitting model (the lowest AICc score) which is referred to as the deltaAIc of 2 or more suggest some support for one model over another and >10 is substantial support. 

logflsmBM<-fitContinuous(tr, logflsm, model="BM")
logflsmWN<-fitContinuous(tr, logflsm, model="white")
logflsmACDC<-fitContinuous(tr, logflsm, model="EB", bounds=list(a=c(-5,5))) # the defaults fix it to only allow an early burst model with decelerating rates (max of the a rate parameter is -0.0000001). So if you want to also fit a late burst model with accelerating rates you have to change the bounds to allow positive numbers.  

logflsmBM$opt$aicc
logflsmWN$opt$aicc
logflsmACDC$opt$aicc

# How can you interpret these results? Think about how the models are related to one another and look at the various parameter estimates by 
logflsmBM$opt
logflsmWN$opt
logflsmACDC$opt

# As clearly noted in the help manual for fitContinuous you need to be aware that the difficulty in finding the optimal solution is determined by an interaction between the nature and complexity of the likelihood space (which is data- and model-dependent) and the numerical optimizer used to explore the space. There is never a guarantee that the optimal solution is found, but using many random starting points (control$niter) and many optimization methods (control$method) will increase these odds."

# Do the conclusions drawn from the contrasts and model-fitting approaches agree?

#==============SECTION 2 part b Changes in BM rate over branches ============== 
# For this section we are going to answer the question of whether dolphins, family Delphinidae have different rates of body size evolution compared to all other cetaceans. First by comparing absolute values of the standardized contrasts within delphinidae (74-95)to all others(1-73 excluding 68 which is the contrast between Delphinidae and other cetaceans)  statistically using the Mann-Whitney U

logmass<-log10(dat[,1])

logmasspic<-pic(logmass, tr, scaled=T, var.contrasts=TRUE)

contrasts4test<-data.frame(abs(logmasspic[,1])[-20]) # this takes the absolute value of the contrasts and removes node 68, which is the contrast between Delphinidae and the other cetaceans. 
contrasts4test$fam<-c(rep(1, 24), rep(2, 22)) # this identifies which contrasts belong to the Delphinidae =2 and the other families = 1, you can identify the relevant nodesusing many different ways but for a fairly small dataset one way is simply to look at the nodelabels (using nodelabels()) after you have plotted the tree) and identify the contrasts for your clade of interest.

boxplot(contrasts4test[,1]~contrasts4test[,2], col=c("darkgrey", "purple"))

wilcox.test(contrasts4test[,1]~contrasts4test[,2]) # this runs the Mann-Whitney U 

# Now we can answer the question by fitting a single-rate of Brownian motion model and comparing it to a two-rate Brownian motion model using the R package OUwie

logmassOUwie<-data.frame(row.names(dat), c(rep(1, 25), rep(2, 23)), log10(dat[,1])) # this sets up a dataset with species rames in the first column, in the second column it identifes the division of the taxa into the different rate categories so in this case it is the identifing all non-delphinid cetaceans as 1 and dolphins as 2, and in the third column we have log10 body mass

logmassBM1<-OUwie(tr, logmassOUwie, model="BM1", simmap.tree=FALSE,  clade=c("Orcinus_orca", "Delphinus_delphis"), diagn=TRUE) #the clade option gives a pair of taxa whose MRCA is the clade of interest, which in this case is Delphinidae, this won't actually matter when model = BM1!
logmassBM1

logmassBMS<-OUwie(tr, logmassOUwie, model=c("BMS"), simmap.tree=FALSE, clade=c("Orcinus_orca", "Delphinus_delphis"), diagn=TRUE) #note root.station=FALSE this stops different different means being estimated for two clades and instead just estimates one at the root
logmassBMS

#What conclusions do you draw and do the results from the PICs and model-fitting agree?

# It is useful to look at the approximate standard errors of sigma^2 (which are calculated from the Hessian matrix)

logmassBMS$solution.se
logmassBM1$solution.se

# Is there are better way of checking the results and the influence of the likelihood surface on the parameter estimates? Yes! it is PARAMETRIC BOOTSTRAPPING! OUwie now has a built in function to run this for you.... this will be covered in later sections. 

