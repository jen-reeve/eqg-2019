require(geiger)   ## Geiger depends on the ape package
data(geospiza)    ## loads dataset geospiza

class(geospiza)
mode(geospiza)
length(geospiza)
names(geospiza)

tree <- geospiza$geospiza.tree  
data <- geospiza$geospiza.data

class(tree)
mode(tree)
length(tree)
names(tree)

plot(tree)
nodelabels()
tiplabels()

tree$edge
class(tree$edge)
mode(tree$edge)
length(tree$edge)
dim(tree$edge)

quartz()   # or x11() or X11() if not on a mac
plot(tree$edge)

?plot.phylo

plot(tree, edge.color="red")
plot(tree, edge.color=c("red", "orange", "yellow"), edge.width=5)
plot(tree, edge.color=c("red", "orange", "yellow"), edge.width=5, edge.lty=2)

tiplabels()
nodelabels()

## Exercise: Plot tree with black branches, make only branch going from 
## node 18 to Platyspiza (11) red. Note that each descendant is only represented 
## once in column 2. This involves making a vector of color names that corresponds 
## to the rows of the edge matrix to supply to edge.color 



## Try coloring or crating a dashed line to whichever branches you choose.








####  make a subtree using the extract.clade function in ape  see ?extract.clade

tree5 <- extract.clade(tree, 21)  ## grab a subtree including node 21 and all descendants
plot(tree5)
nodelabels()
tiplabels()

bl5 <-  tree5$edge.length   ## save the original branch lengths as bl5
bl <- rnorm(length(bl5), mean=mean(bl5), sd=sd(bl5))   ## random branch lengths

### Exercise: replace the branch lengths in your subtree with the random ones




### Exercise: create functions and use them to do some simulations of random branch lengths

bl <- -1
while(any(bl < 0)) bl <- rnorm(length(bl5), mean=mean(bl5), sd=sd(bl5))

gen.bl <- function( ttree ) {
	bl <- -1
	el <- ttree$edge.length
	while(any(bl < 0)) bl <- rnorm(length(el), mean=mean(el), sd=sd(el))
	return(bl)
}

change.bl <- function( ttree, bl ){
	ttree$edge.length <- bl
	return(ttree)
}

branchlengths <- gen.bl(tree5)
simtree <- change.bl(tree5, branchlengths)
plot(simtree)


### Exercise: Branch length simulation 2. Instead of drawing each branch length from a random normal distribution, letʻs assume that the original branch lengths are estimated with some error, letʻs assume that the sd = 25% of the branch length value. 

	