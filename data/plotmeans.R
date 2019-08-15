
medianize <- function(ev) {
# "medianize" principal component (ev) by getting the medians
# of the values for the  x  coordinates or for the  y  coordinates
p <- length(ev);
medx <- median(ev[seq(1,(p-1),2)]);
medy <- median(ev[seq(2,p,2)]);
medxy <- medx*rep(c(1,0),p/2)+medy*rep(c(0,1),p/2);
return(medxy);
}


plotmeansandpc <- function(a, pc) {
# plot mean form, with PC loadings as red arrows
  n1 <- dim(a)[1];                # number of forms in array a
  p1 <- dim(a)[2]-1;              # number of coordinates (not counting name)
  if ((dim(a)[2] %% 2) == 0) {
    p1 <- p1 - 1;                 # ... and not counting log-size if it is there
  }
  minx <- min(a[ ,seq(2,p1,2)]);      # minimum of x's
  maxx <- max(a[ ,seq(2,p1,2)]);      # maximum of x's
  miny <- min(a[ ,seq(3,(p1+1),2)]);    # minimum of y's
  maxy <- max(a[ ,seq(3,(p1+1),2)]);    # maximum of y's
  btemp <- colSums(a[,2:(p1+1)])/n1;    # get coordinates of the mean form
  b <- matrix(c(btemp[1:p1],btemp[1],btemp[2]),2,(p1/2)+1);
                                           #  coordinates for outline of mean form
                                           #  with starting point duplicated at end
  plot(t(b), type = "l", asp = 1, xlab = "", ylab = "", xlim = c(minx, maxx), ylim = c(miny, maxy), col="gray",lwd=4);
                                           # plot mean form as gray outline
  ctemp <- matrix(pc, 2, p1/2);            # get a matrix of the PC
  for (i in 1:(p1/2)) {                    # plot red arrows for pc on mean form landmarks
    arrows( as.numeric(btemp[2*i-1]) - as.numeric(ctemp[1, i]), 
            as.numeric(btemp[2*i]) - as.numeric(ctemp[2, i]),
            as.numeric(btemp[2*i-1]) + as.numeric(ctemp[1, i]), 
            as.numeric(btemp[2*i]) + as.numeric(ctemp[2, i]),
            length = 0.05, angle = 20, col = "red", lwd = 4);
  }
}


plotmeansandpcmed <- function(z,pc) {
# plot mean form in gray and PCs in red for departure from medianized PC
  plotmeansandpc(z, pc-medianize(pc));
}



