# Simulation 

# There are two interelated questions. The first is given the tree and the data do we have the power to distinguish between the different models we are fitting? The second is what is the uncertainty associated with the estimated parameters for the models? This exercise deals with the latter.

# We simulate under our models with the estimated parameter values and re-estimate the model-fit and parameters. This can tell us if we have the power to distinguish between the models we are fitting. If we do have the power we should be able to simulate data under the estimated parameters for each of the models and to be able to distinguish between them. If we don't then we make inaccurate biological statements based on estimates from analyses that have inadequate power. It can also tell us the uncertainty associated with the parameter estimates.


## 2) Using the integrated OUwie bootstrap function to calculate 95% Confidence Intervals around the parameter estimates from the best-fitting model 
# This will provide an estimate of the degree of uncertainty in our parameter estimates including the influence that the likelihood surface has on the parameter estimates and is a much more rigorous way to present the results compared to the single value estimated from each phylogeny. We are testing the hypothesis that diet influences body mass evolution across marsupials, the expectation is that herbivores (1) will have larger optima than omnivores (2) or carnivores (3).

library(OUwie)
dietmtr<-read.nexus("marsupialdiettree.nex")
marsupiald<-read.table("marsupialmass_diet.txt", header=T)

mods<-c("BM1", "OU1", "OUM", "OUMA", "OUMV") # pick all the models you want to fit - we are using all the OU models that allow different for diet as that matches our prediction and also the two simpler models

sizedietresults<-list()# setting up an empty list for the results
for(i in 1:length(mods)){#loop to run OUwie across all the models you want to fit
	sizedietresults[[i]]<-OUwie(dietmtr, marsupiald, model=mods[i], simmap.tree=FALSE, root.station=TRUE, diagn=T) #note you probably don't actually want to run this with root.station=TRUE for the BMS model
}
dietbestfit<-c(sizedietresults[[1]]$AICc, sizedietresults[[2]]$AICc, sizedietresults[[3]]$AICc, sizedietresults[[4]]$AICc, sizedietresults[[5]]$AICc)
names(dietbestfit)<-mods# 
dietbestfit-min(dietbestfit)

# We want to use the parameters from the best-fitting model
sizedietresults[[4]]

# This will take a long time to run but you can read in the results from the analyses of 1000 I ran earlier (takes about 8 hrs on my mac 3.1 GHz Intel Core i7)

sizedietoumboot1000<-OUwie.boot(dietmtr, marsupiald, model=c("OUMA"), nboot=1000, alpha=sizedietresults[[4]]$solution[1,] , sigma.sq=sizedietresults[[4]]$solution[2,], theta=sizedietresults[[4]]$theta[,1], theta0= sizedietresults[[4]]$theta[2,1]) #When theta0 was dropped from the model (theta0=T) i.e. not estimated, then we need to set the root as the the value of the continuous trait for the selective regime mapped at the root which is omnivore (2)

write.table(sizedietoumboot1000, file="sizedietoumboot1000.txt")

sizedietoumboot1000 <-read.table("sizedietoumboot1000.txt")

head(sizedietoumboot1000)

# Plot the entire distribution of the alpha parameter (sigma was fixed between the different diets as the best fitting model was OUMA)
plot(density(sizedietoumboot1000[,1]), col="green", ylim=c(0,30), xlab="Marsupial Size Alpha Estimate") # distribution of the alpha estimates for herbivores
lines(density(sizedietoumboot1000[,2]), col="purple")# distribution of the alpha estimates for omnivores
lines(density(sizedietoumboot1000[,3]), col="blue")# distribution of the alpha estimates for carnivores
abline(v= sizedietresults[[4]]$solution[1,1], col="green")# This is the empirical estimate of alpha
abline(v= sizedietresults[[4]]$solution[1,2], col="purple")# This is the empirical estimate of alpha
abline(v= sizedietresults[[4]]$solution[1,3], col="blue")# This is the empirical estimate of alpha

# We can also calculate the 95% CI of the parameter estimates easily.
C195sizedietoumboot1000 <-apply(sizedietoumboot1000, 2, quantile, probs=c(0.025,0.975)) 

# What conclusions can you draw from this? Is it different for the theta parameter?

# Why might you not use bootstrapping?  If your dataset is very large then the time taken to run the simulations may be prohibitive, especially if you are running it on 1000 stochastically mapped trees.


## 4)Estimate the 95% CI interval around the estimate of Pagel's lambda for parrotfish. Pagel's lambda is frequently used to estimate the degree of phylogenetic signal within a continuous trait or of the residuals from a phylogenetic regression by estimating the optimal transformation of the internal branches. Thus lambda=0 you have a star phylogeny and no signal but when lambda=1 is the identical phylogeny and indicative of Brownian motion. Pagel's lambda is estimated in geiger's fitContinuous function model="lambda"


library(geiger)

ptr<-read.nexus("parrotfishtree.nex")
parrotfishmorph<-read.table("parrotfishesRawMorphology.txt", header=T, row.names=1)
rownames(parrotfishmorph)[1]<-"Bolbometopon_muricatum" #this is a misspelling in the dataset!

foo<-treedata(ptr, parrotfishmorph, sort=T) # match dataset and tree

prunedmorph<-foo[[2]]
prunedtr<-foo[[1]]

loglength<-log(prunedmorph[,1])
names(loglength)<-rownames(prunedmorph)

lambdalength<-fitContinuous(prunedtr, loglength, model="lambda") #Estimate the ML estimate of lambda, root state z0 and sigma^2

lambdatree<-rescale(prunedtr, model="lambda",  lambdalength$opt$lambda) # to simulate under the lambda model we transform the branch lengths of the tree by the estimated lambda parameter and then run a simple Brownian motion on the transformed tree

lambdasim<-sim.char(lambdatree, lambdalength$opt$sigsq, nsim=500, model="BM", root= lambdalength$opt$z0)

lambdasimres<-vector()
for(i in 1:dim(lambdasim)[3]){ # a loop that fits the lambda model to the 500 simulated characters using fitContinuous and puts the estimated lambda value in a vector called lambdasimres
	print(i)
	tmp<-fitContinuous(prunedtr , lambdasim[,,i], model="lambda")
	lambdasimres<-c(lambdasimres, tmp$opt$lambda)
}

plot(density(lambdasimres))
abline(v= lambdalength$opt$lambda)
parrotfishlengthlambdaCI<-quantile(lambdasimres, c(0.025, 0.975))

## What conclusions do you draw about the uncertainty when estimating lambda for this dataset? What does this tell us about the phylogenetic signal within the data?